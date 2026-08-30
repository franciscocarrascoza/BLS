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

bool Analyzer::processFrame(const Frame& frame, FrameMetrics& metrics, std::string& err,
                            std::vector<int>* labels) {
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

  // Helper lambda to compute coordinate bounds
  auto computeBounds = [&](Vec3& minPos, Vec3& maxPos) {
    minPos = Vec3{std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity()};
    maxPos = Vec3{-std::numeric_limits<double>::infinity(),
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
  };

  // Check if box needs correction (zero/invalid or unreasonably large)
  bool needsBoxCorrection = (len0 < 1e-8 || len1 < 1e-8 || len2 < 1e-8);

  if (!needsBoxCorrection && !frame.xyz.empty() && config_.boxMode == BoxMode::Auto) {
    // Box appears valid, but check if it's unreasonably large compared to coordinates
    Vec3 minPos, maxPos;
    computeBounds(minPos, maxPos);

    double coordExtentX = maxPos.x - minPos.x;
    double coordExtentY = maxPos.y - minPos.y;
    double coordExtentZ = maxPos.z - minPos.z;

    // If box is more than 10x larger than coordinate extent in any dimension, it's likely incorrect
    const double suspiciousRatio = 10.0;
    bool boxTooLarge = (len0 > coordExtentX * suspiciousRatio) ||
                       (len1 > coordExtentY * suspiciousRatio) ||
                       (len2 > coordExtentZ * suspiciousRatio);

    if (boxTooLarge) {
      Logger::warn("Box size (", len0, " x ", len1, " x ", len2,
                   ") is unreasonably large compared to coordinate extent (",
                   coordExtentX, " x ", coordExtentY, " x ", coordExtentZ,
                   "). Auto-correcting to fit coordinates.");
      needsBoxCorrection = true;
    }
  }

  if (needsBoxCorrection && !frame.xyz.empty()) {
    Vec3 minPos, maxPos;
    computeBounds(minPos, maxPos);
    const double padding = config_.gridSpacing * 2.0;
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

  // Check memory requirements before allocation
  std::size_t requiredMemory = estimateGridMemoryBytes(nx, ny, nz);
  std::size_t availableMemory = availableSystemRAMBytes();

  // Use 80% of available RAM as the safety limit
  std::size_t maxAllowedMemory = static_cast<std::size_t>(availableMemory * 0.8);

  if (requiredMemory > maxAllowedMemory) {
    std::ostringstream oss;
    oss << "Grid allocation would require " << (requiredMemory / (1024.0 * 1024.0 * 1024.0))
        << " GB, which exceeds available RAM limit ("
        << (maxAllowedMemory / (1024.0 * 1024.0 * 1024.0)) << " GB).\n";
    oss << "Grid dimensions: " << nx << " x " << ny << " x " << nz << " = "
        << (static_cast<std::size_t>(nx) * ny * nz) << " voxels\n";
    oss << "Box size: " << norm(col0) << " x " << norm(col1) << " x " << norm(col2) << " Angstroms\n";
    oss << "Grid spacing: " << config_.gridSpacing << " Angstroms\n\n";
    oss << "Solutions:\n";
    oss << "  1. Increase GRID_SPACING (current: " << config_.gridSpacing << " A)\n";
    double minSpacing = std::max({norm(col0), norm(col1), norm(col2)}) /
                        maxBoxDimensionForRAM(1.0, maxAllowedMemory);
    oss << "     Minimum spacing for this box: " << minSpacing << " A\n";
    oss << "  2. Reduce box size (check CRYST1 record in PDB or use BOX MANUAL in config)\n";
    double maxBoxSize = maxBoxDimensionForRAM(config_.gridSpacing, maxAllowedMemory);
    oss << "     Maximum box dimension for current spacing: " << maxBoxSize << " A";
    err = oss.str();
    return false;
  }

  if (!impl_->configured || nx != impl_->grid.nx() || ny != impl_->grid.ny() ||
      nz != impl_->grid.nz()) {
    impl_->grid.configure(nx, ny, nz, config_.gridSpacing, activeBox, origin,
                          periodicityForBoxMode(config_.boxMode));
    impl_->configured = true;
  } else {
    impl_->grid.reset();
  }

  impl_->grid.rasterize(frame.xyz, selectionPtr, config_.cutoff, config_.occupancy);

  Enumerator enumerator(impl_->scaledBasis, impl_->lattice.offsets, nx, ny, nz,
                        impl_->grid.occupancy());

  SkipDFSConfig skipCfg{nx, ny, nz, config_.connectivity, config_.refinementStride};
  SkipDFS dfs(skipCfg, impl_->grid.occupancy(), impl_->grid.visited());

  int seeds = 0;
  int seedHits = 0;
  int nclusters = 0;
  int maxCluster = 0;
  std::size_t refinedVoxels = 0;
  std::vector<int> clusterSizes;

  const std::vector<uint8_t>& occ = impl_->grid.occupancy();
  std::vector<uint8_t>& visited = impl_->grid.visited();
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
