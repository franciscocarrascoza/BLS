#include "bls/BLS.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "grid/Grid.hpp"
#include "grid/GridSpec.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"
#include "refine/SkipDFS.hpp"
#include "refine/StandardDFS.hpp"
#include "util/Logging.hpp"
#include "util/RSS.hpp"
#include "util/Timer.hpp"

namespace bls {

namespace {

double computeDnnVoxel(const BLSConfig& config) {
  if (config.hasExplicitDnn && config.dnn > 0.0) {
    return config.dnn;
  }

  if (!config.radii.empty()) {
    double minSum = std::numeric_limits<double>::infinity();
    for (double ri : config.radii) {
      for (double rj : config.radii) {
        minSum = std::min(minSum, ri + rj);
      }
    }
    if (std::isfinite(minSum) && minSum > 0.0) {
      double minSumVoxel = minSum / config.gridSpacing;
      return config.alpha * minSumVoxel;
    }
  }

  return std::ceil(1.0 / config.gridSpacing);
}

}  // namespace

struct Analyzer::Impl {
  explicit Impl(const BLSConfig& cfg)
      : lattice(buildLattice(cfg.lattice)),
        dnnVoxel(computeDnnVoxel(cfg)),
        scaledBasis(lattice.basis * (dnnVoxel / lattice.dmin)),
        latticeSettings(cfg.lattice),
        latticeName(latticeToString(cfg.lattice.lattice)),
        centeringName(centeringToString(cfg.lattice.centering)) {}

  Grid grid;
  LatticeDescriptor lattice;
  double dnnVoxel;
  Mat3 scaledBasis;
  LatticeSettings latticeSettings;
  std::string latticeName;
  std::string centeringName;
};

Analyzer::Analyzer(const BLSConfig& config) : config_(config), impl_(new Impl(config)) {}

Analyzer::~Analyzer() = default;

void Analyzer::setSelection(const std::vector<int>& indices, int natoms) {
  selection_ = indices;
  selectionIsAll_ = selection_.empty();
  if (!selectionIsAll_) {
    for (int idx : selection_) {
      if (idx < 0 || idx >= natoms) {
        throw std::runtime_error("Selection index outside topology atom count.");
      }
    }
  }
}

bool Analyzer::processFrame(const Frame& frame, FrameMetrics& metrics, std::string& err,
                            std::vector<int>* labels) {
  ScopedTimer totalTimer;
  const std::vector<int>* selectionPtr = selectionIsAll_ ? nullptr : &selection_;
  GridSpec spec;
  if (!deriveGridSpec(config_, frame, selectionPtr, spec, err)) return false;
  configureGrid(impl_->grid, spec);
  impl_->grid.rasterize(frame.xyz, selectionPtr, config_.cutoff, config_.occupancy);
  if (!labelGrid(impl_->grid, metrics, err, labels)) return false;
  metrics.totalMs = totalTimer.elapsedMilliseconds();
  metrics.timePs = frame.time;
  metrics.natoms = frame.natoms;
  metrics.nx = spec.nx;
  metrics.ny = spec.ny;
  metrics.nz = spec.nz;
  return true;
}

bool Analyzer::probeBasis(int nx, int ny, int nz, bool periodic, Mat3& basis, double& dnnVoxel,
                          std::string& err) const {
  basis = impl_->scaledBasis;
  dnnVoxel = impl_->dnnVoxel;
  if (!periodic) return true;
  return commensurateCubicBasis(impl_->latticeSettings, impl_->lattice, impl_->dnnVoxel, nx, ny,
                                nz, basis, dnnVoxel, err);
}

bool Analyzer::labelGrid(Grid& grid, FrameMetrics& metrics, std::string& err,
                         std::vector<int>* labels) {
  ScopedTimer timer;
  const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
  const std::vector<uint8_t>& occ = grid.occupancy();
  std::vector<uint8_t>& visited = grid.visited();
  // BLS's refinement needs a zeroed visited array, as DFS does; it is initialised here,
  // inside the labelling timer, exactly as every comparison method initialises its own
  // (e.g. traditionalDFS, Algorithms.cpp). Relying on the harness's zeroing outside the
  // timer would give BLS alone a free grid-sized clear (audit S1-2).
  std::fill(visited.begin(), visited.end(), 0);
  const double clearEnd = timer.elapsedMilliseconds();

  // Probe evaluation (§2.1.2, audit D3): the lattice sites inside the grid are evaluated
  // and the occupied ones seed the refinement. Under PBC the lattice is first made
  // commensurate with the periodic grid (audit D10) -- per frame, since the grid
  // dimensions can change between frames -- and the reported dNN is the effective one.
  const bool periodic = grid.periodicity() == BoxPeriodicity::Periodic;
  double dnnVoxel = 0.0;
  Mat3 basis;
  if (!probeBasis(nx, ny, nz, periodic, basis, dnnVoxel, err)) return false;
  const Vec3 latticeOrigin{config_.latticeOrigin[0], config_.latticeOrigin[1],
                           config_.latticeOrigin[2]};  // deck LATTICE_ORIGIN, voxels
  const double probeStart = timer.elapsedMilliseconds();
  const Enumerator enumerator =
      periodic ? Enumerator(basis, impl_->lattice.offsets, nx, ny, nz, occ,
                            Enumerator::PeriodicTag{}, latticeOrigin)
               : Enumerator(basis, impl_->lattice.offsets, nx, ny, nz, occ, latticeOrigin);
  const double probeEnd = timer.elapsedMilliseconds();

  int seeds = 0;
  int seedHits = 0;
  int nclusters = 0;
  int maxCluster = 0;
  std::size_t refinedVoxels = 0;
  std::vector<int> clusterSizes;

  if (labels) {
    labels->assign(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                       static_cast<std::size_t>(nz),
                   -1);
  }

  // Refinement (stage 2) from every seed. The seeds are occupied by construction (the
  // Enumerator keeps only occupied voxels), and runFrom returns 0 for a voxel an earlier walk
  // has already claimed, so no check is repeated here.
  const SkipDFSConfig refineCfg{nx, ny, nz, config_.connectivity, config_.refinementStride,
                                periodic};
  auto refineFromSeeds = [&](auto& refiner) {
    enumerator.forEach([&](const Enumerator::Seed& seed) {
      ++seeds;
      // nclusters is the count of components already accepted, so it is the
      // 0-based ordinal of this one: dense by construction.
      const int size = refiner.runFrom(seed.x, seed.y, seed.z, labels, nclusters);
      if (size > 0) {
        ++seedHits;
        ++nclusters;
        maxCluster = std::max(maxCluster, size);
        clusterSizes.push_back(size);
        refinedVoxels += refiner.refinedVoxels();
      }
    });
  };
  if (config_.refinement == RefinementMode::DFS) {
    StandardDFS refiner(refineCfg, occ, visited);
    refineFromSeeds(refiner);
  } else {
    SkipDFS refiner(refineCfg, occ, visited);
    refineFromSeeds(refiner);
  }
  const double refineEnd = timer.elapsedMilliseconds();

  metrics.dnnVoxel = dnnVoxel;
  metrics.lattice = impl_->latticeName;
  metrics.centering = impl_->centeringName;
  metrics.seeds = seeds;
  metrics.seedHits = seedHits;
  metrics.nclusters = nclusters;
  metrics.maxCluster = maxCluster;
  metrics.refinedVoxels = refinedVoxels;
  metrics.probes = enumerator.probes();
  std::sort(clusterSizes.begin(), clusterSizes.end(), std::greater<int>());
  metrics.clusterSizes = std::move(clusterSizes);
  metrics.elapsedMs = timer.elapsedMilliseconds();
  metrics.clearMs = clearEnd;
  metrics.probeMs = probeEnd - probeStart;
  metrics.refineMs = refineEnd - probeEnd;
  metrics.refinement = config_.refinement == RefinementMode::DFS ? "dfs" : "skip_dfs";
  metrics.nx = nx;
  metrics.ny = ny;
  metrics.nz = nz;
  return true;
}

}  // namespace bls
