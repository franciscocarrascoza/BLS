#include "bls/BLS.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "grid/Grid.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"
#include "refine/SkipDFS.hpp"
#include "util/Logging.hpp"
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

Mat3 makeManualBox(const ManualBox& box) {
  double lx = box.xhi - box.xlo;
  double ly = box.yhi - box.ylo;
  double lz = box.zhi - box.zlo;
  if (lx <= 0 || ly <= 0 || lz <= 0) {
    throw std::runtime_error("Manual box extents must be positive.");
  }
  return Mat3{Vec3{lx, 0.0, 0.0}, Vec3{0.0, ly, 0.0}, Vec3{0.0, 0.0, lz}};
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
  bool configured{false};
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

bool Analyzer::processFrame(const Frame& frame, FrameMetrics& metrics, std::string& err) {
  ScopedTimer timer;

  Mat3 activeBox;
  Vec3 origin{0.0, 0.0, 0.0};
  if (config_.boxMode == BoxMode::Manual) {
    try {
      activeBox = makeManualBox(config_.manualBox);
    } catch (const std::exception& ex) {
      err = ex.what();
      return false;
    }
    origin = Vec3{config_.manualBox.xlo, config_.manualBox.ylo, config_.manualBox.zlo};
  } else {
    activeBox = frame.box;
  }

  auto col0 = activeBox.column(0);
  auto col1 = activeBox.column(1);
  auto col2 = activeBox.column(2);

  const std::vector<int>* selectionPtr = selectionIsAll_ ? nullptr : &selection_;

  const double len0 = norm(col0);
  const double len1 = norm(col1);
  const double len2 = norm(col2);
  if ((len0 < 1e-8 || len1 < 1e-8 || len2 < 1e-8) && !frame.xyz.empty()) {
    Vec3 minPos{std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    Vec3 maxPos{-std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity()};
    auto accumulate = [&](const Vec3& p) {
      minPos.x = std::min(minPos.x, p.x);
      minPos.y = std::min(minPos.y, p.y);
      minPos.z = std::min(minPos.z, p.z);
      maxPos.x = std::max(maxPos.x, p.x);
      maxPos.y = std::max(maxPos.y, p.y);
      maxPos.z = std::max(maxPos.z, p.z);
    };
    if (selectionPtr) {
      for (int idx : *selectionPtr) {
        if (idx >= 0 && idx < frame.natoms) accumulate(frame.xyz[static_cast<std::size_t>(idx)]);
      }
    } else {
      for (const auto& p : frame.xyz) accumulate(p);
    }
    const double padding = config_.gridSpacing;
    origin = minPos;
    activeBox = Mat3{Vec3{std::max(maxPos.x - minPos.x + padding, padding), 0.0, 0.0},
                     Vec3{0.0, std::max(maxPos.y - minPos.y + padding, padding), 0.0},
                     Vec3{0.0, 0.0, std::max(maxPos.z - minPos.z + padding, padding)}};
    col0 = activeBox.column(0);
    col1 = activeBox.column(1);
    col2 = activeBox.column(2);
  }
  int nx = std::max(1, static_cast<int>(std::ceil(norm(col0) / config_.gridSpacing)));
  int ny = std::max(1, static_cast<int>(std::ceil(norm(col1) / config_.gridSpacing)));
  int nz = std::max(1, static_cast<int>(std::ceil(norm(col2) / config_.gridSpacing)));

  if (!impl_->configured || nx != impl_->grid.nx() || ny != impl_->grid.ny() ||
      nz != impl_->grid.nz()) {
    impl_->grid.configure(nx, ny, nz, config_.gridSpacing, activeBox, origin);
    impl_->configured = true;
  } else {
    impl_->grid.reset();
  }

  impl_->grid.rasterize(frame.xyz, selectionPtr, config_.cutoff, config_.occupancy);

  Enumerator enumerator(impl_->scaledBasis, impl_->lattice.offsets, nx, ny, nz);

  SkipDFSConfig skipCfg{nx, ny, nz, config_.connectivity, config_.skip};
  SkipDFS dfs(skipCfg, impl_->grid.occupancy(), impl_->grid.visited());

  int seeds = 0;
  int seedHits = 0;
  int nclusters = 0;
  int maxCluster = 0;
  std::size_t refinedVoxels = 0;
  std::vector<int> clusterSizes;

  const std::vector<uint8_t>& occ = impl_->grid.occupancy();
  std::vector<uint8_t>& visited = impl_->grid.visited();

  enumerator.forEach([&](const Enumerator::Seed& seed) {
    ++seeds;
    std::size_t idx =
        static_cast<std::size_t>(seed.x) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz) +
        static_cast<std::size_t>(seed.y) * static_cast<std::size_t>(nz) +
        static_cast<std::size_t>(seed.z);
    if (!occ[idx] || visited[idx]) {
      return;
    }
    int size = dfs.runFrom(seed.x, seed.y, seed.z);
    if (size > 0) {
      ++seedHits;
      ++nclusters;
      maxCluster = std::max(maxCluster, size);
      clusterSizes.push_back(size);
      refinedVoxels += dfs.refinedVoxels();
    }
  });

  metrics.timePs = frame.time;
  metrics.natoms = frame.natoms;
  metrics.nx = nx;
  metrics.ny = ny;
  metrics.nz = nz;
  metrics.dnnVoxel = impl_->dnnVoxel;
  metrics.lattice = impl_->latticeName;
  metrics.centering = impl_->centeringName;
  metrics.seeds = seeds;
  metrics.seedHits = seedHits;
  metrics.nclusters = nclusters;
  metrics.maxCluster = maxCluster;
  metrics.refinedVoxels = refinedVoxels;
  std::sort(clusterSizes.begin(), clusterSizes.end(), std::greater<int>());
  metrics.clusterSizes = std::move(clusterSizes);
  metrics.elapsedMs = timer.elapsedMilliseconds();

  return true;
}

}  // namespace bls
