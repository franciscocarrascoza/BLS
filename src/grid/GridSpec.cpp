// Created: 2026-09-26T16:40+02:00 | by: CC audit | purpose: the one shared grid builder (fixes Defect 12)
#include "grid/GridSpec.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "grid/Grid.hpp"
#include "util/Logging.hpp"
#include "util/RSS.hpp"

namespace bls {

double GridSpec::h(int axis) const {
  const int n = axis == 0 ? nx : (axis == 1 ? ny : nz);
  return norm(box.column(axis)) / static_cast<double>(n);
}

namespace {

void selectedBounds(const Frame& frame, const std::vector<int>* selection, Vec3& minPos,
                    Vec3& maxPos) {
  minPos = Vec3{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
  maxPos = Vec3{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity()};
  auto accumulate = [&](const Vec3& p) {
    minPos.x = std::min(minPos.x, p.x);
    minPos.y = std::min(minPos.y, p.y);
    minPos.z = std::min(minPos.z, p.z);
    maxPos.x = std::max(maxPos.x, p.x);
    maxPos.y = std::max(maxPos.y, p.y);
    maxPos.z = std::max(maxPos.z, p.z);
  };
  if (selection) {
    for (int idx : *selection) {
      if (idx >= 0 && idx < frame.natoms) accumulate(frame.xyz[static_cast<std::size_t>(idx)]);
    }
  } else {
    for (const auto& p : frame.xyz) accumulate(p);
  }
}

}  // namespace

bool deriveGridSpec(const BLSConfig& config, const Frame& frame,
                    const std::vector<int>* selection, GridSpec& spec, std::string& err) {
  spec = GridSpec{};
  spec.spacing = config.gridSpacing;
  spec.pbc = config.pbc;

  Mat3 activeBox = frame.box;
  Vec3 origin{0.0, 0.0, 0.0};

  if (config.boxMode == BoxMode::Manual) {
    const ManualBox& b = config.manualBox;
    const double lx = b.xhi - b.xlo, ly = b.yhi - b.ylo, lz = b.zhi - b.zlo;
    if (lx <= 0 || ly <= 0 || lz <= 0) {
      err = "Manual box extents must be positive.";
      return false;
    }
    activeBox = Mat3{Vec3{lx, 0.0, 0.0}, Vec3{0.0, ly, 0.0}, Vec3{0.0, 0.0, lz}};
    origin = Vec3{b.xlo, b.ylo, b.zlo};
  } else if (config.boxMode == BoxMode::Cell) {
    if (norm(activeBox.column(0)) < 1e-8 || norm(activeBox.column(1)) < 1e-8 ||
        norm(activeBox.column(2)) < 1e-8 || std::fabs(determinant(activeBox)) < 1e-8) {
      err = "BOX CELL: the frame carries no usable cell (no CRYST1 record, or a degenerate one).";
      return false;
    }
  } else {
    // BOX AUTO -- identical arithmetic to the pre-audit Analyzer::processFrame.
    const double len0 = norm(activeBox.column(0));
    const double len1 = norm(activeBox.column(1));
    const double len2 = norm(activeBox.column(2));
    bool needsBoxCorrection = (len0 < 1e-8 || len1 < 1e-8 || len2 < 1e-8);
    if (!needsBoxCorrection && !frame.xyz.empty()) {
      Vec3 minPos, maxPos;
      selectedBounds(frame, selection, minPos, maxPos);
      const double ex = maxPos.x - minPos.x, ey = maxPos.y - minPos.y, ez = maxPos.z - minPos.z;
      const double suspiciousRatio = 10.0;
      if (len0 > ex * suspiciousRatio || len1 > ey * suspiciousRatio ||
          len2 > ez * suspiciousRatio) {
        Logger::warn("Box size (", len0, " x ", len1, " x ", len2,
                     ") is unreasonably large compared to coordinate extent (", ex, " x ", ey,
                     " x ", ez, "). Auto-correcting to fit coordinates.");
        needsBoxCorrection = true;
      }
    }
    if (needsBoxCorrection && !frame.xyz.empty()) {
      Vec3 minPos, maxPos;
      selectedBounds(frame, selection, minPos, maxPos);
      const double padding = config.gridSpacing * 2.0;
      origin = minPos;
      activeBox = Mat3{Vec3{std::max(maxPos.x - minPos.x + padding, padding), 0.0, 0.0},
                       Vec3{0.0, std::max(maxPos.y - minPos.y + padding, padding), 0.0},
                       Vec3{0.0, 0.0, std::max(maxPos.z - minPos.z + padding, padding)}};
    }
  }

  const Vec3 col0 = activeBox.column(0), col1 = activeBox.column(1), col2 = activeBox.column(2);
  const int nx = std::max(1, static_cast<int>(std::ceil(norm(col0) / config.gridSpacing)));
  const int ny = std::max(1, static_cast<int>(std::ceil(norm(col1) / config.gridSpacing)));
  const int nz = std::max(1, static_cast<int>(std::ceil(norm(col2) / config.gridSpacing)));

  const std::size_t requiredMemory = estimateGridMemoryBytes(nx, ny, nz);
  const std::size_t maxAllowedMemory = static_cast<std::size_t>(availableSystemRAMBytes() * 0.8);
  if (requiredMemory > maxAllowedMemory) {
    std::ostringstream oss;
    oss << "Grid allocation would require " << (requiredMemory / (1024.0 * 1024.0 * 1024.0))
        << " GB, which exceeds available RAM limit ("
        << (maxAllowedMemory / (1024.0 * 1024.0 * 1024.0)) << " GB).\n";
    oss << "Grid dimensions: " << nx << " x " << ny << " x " << nz << " = "
        << (static_cast<std::size_t>(nx) * ny * nz) << " voxels\n";
    oss << "Box size: " << norm(col0) << " x " << norm(col1) << " x " << norm(col2)
        << " Angstroms\n";
    oss << "Grid spacing: " << config.gridSpacing << " Angstroms\n\n";
    oss << "Solutions:\n";
    oss << "  1. Increase GRID_SPACING (current: " << config.gridSpacing << " A)\n";
    const double minSpacing = std::max({norm(col0), norm(col1), norm(col2)}) /
                              maxBoxDimensionForRAM(1.0, maxAllowedMemory);
    oss << "     Minimum spacing for this box: " << minSpacing << " A\n";
    oss << "  2. Reduce box size (check CRYST1 record in PDB or use BOX MANUAL in config)\n";
    oss << "     Maximum box dimension for current spacing: "
        << maxBoxDimensionForRAM(config.gridSpacing, maxAllowedMemory) << " A";
    err = oss.str();
    return false;
  }

  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.box = activeBox;
  spec.origin = origin;
  return true;
}

void configureGrid(Grid& grid, const GridSpec& spec) {
  grid.configure(spec.nx, spec.ny, spec.nz, spec.spacing, spec.box, spec.origin,
                 spec.periodicity());
}

}  // namespace bls
