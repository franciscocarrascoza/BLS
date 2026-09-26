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
        latticeName(latticeToString(cfg.lattice.lattice)),
        centeringName(centeringToString(cfg.lattice.centering)) {}

  Grid grid;
  LatticeDescriptor lattice;
  double dnnVoxel;
  Mat3 scaledBasis;
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

bool Analyzer::labelGrid(Grid& grid, FrameMetrics& metrics, std::string& err,
                         std::vector<int>* labels) {
  (void)err;
  ScopedTimer timer;
  const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
  const std::vector<uint8_t>& occ = grid.occupancy();
  std::vector<uint8_t>& visited = grid.visited();
  // BLS's refinement needs a zeroed visited array, as DFS does; it is initialised here,
  // inside the labelling timer, exactly as every comparison method initialises its own
  // (e.g. traditionalDFS, Algorithms.cpp). Relying on the harness's zeroing outside the
  // timer would give BLS alone a free grid-sized clear (audit S1-2).
  std::fill(visited.begin(), visited.end(), 0);

  Enumerator enumerator(impl_->scaledBasis, impl_->lattice.offsets, nx, ny, nz, occ);

  SkipDFSConfig skipCfg{nx, ny, nz, config_.connectivity, config_.refinementStride};
  SkipDFS dfs(skipCfg, occ, visited);

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

  enumerator.forEach([&](const Enumerator::Seed& seed) {
    ++seeds;
    std::size_t idx =
        static_cast<std::size_t>(seed.x) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz) +
        static_cast<std::size_t>(seed.y) * static_cast<std::size_t>(nz) +
        static_cast<std::size_t>(seed.z);
    if (!occ[idx] || visited[idx]) {
      return;
    }
    // nclusters is the count of components already accepted, so it is the
    // 0-based ordinal of this one: dense by construction.
    int size = dfs.runFrom(seed.x, seed.y, seed.z, labels, nclusters);
    if (size > 0) {
      ++seedHits;
      ++nclusters;
      maxCluster = std::max(maxCluster, size);
      clusterSizes.push_back(size);
      refinedVoxels += dfs.refinedVoxels();
    }
  });

  metrics.dnnVoxel = impl_->dnnVoxel;
  metrics.lattice = impl_->latticeName;
  metrics.centering = impl_->centeringName;
  metrics.seeds = seeds;
  metrics.seedHits = seedHits;
  metrics.nclusters = nclusters;
  metrics.maxCluster = maxCluster;
  metrics.refinedVoxels = refinedVoxels;
  metrics.probes = 0;
  std::sort(clusterSizes.begin(), clusterSizes.end(), std::greater<int>());
  metrics.clusterSizes = std::move(clusterSizes);
  metrics.elapsedMs = timer.elapsedMilliseconds();
  metrics.nx = nx;
  metrics.ny = ny;
  metrics.nz = nz;
  return true;
}

}  // namespace bls
