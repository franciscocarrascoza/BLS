// Created: 2026-09-26T16:40+02:00 | by: CC audit | purpose: the one shared grid builder (fixes Defect 12: BLS and the comparison methods derived their boxes separately and main.cpp never read BOX MANUAL)
#pragma once

#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "common/Types.hpp"
#include "io/TrajectoryReader.hpp"

namespace bls {

class Grid;

// Everything that defines the grid a frame is labelled on. Derived once per frame by
// deriveGridSpec() and consumed unchanged by BLS and by every comparison method; no
// method re-derives box, origin, dimensions or periodicity.
struct GridSpec {
  int nx{0}, ny{0}, nz{0};
  double spacing{0.0};  // nominal GRID_SPACING (Angstrom)
  Mat3 box{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};  // columns = cell vectors (Angstrom)
  Vec3 origin{0.0, 0.0, 0.0};                              // Angstrom
  PeriodicAxes pbc{};

  BoxPeriodicity periodicity() const { return periodicityFor(pbc); }
  // Actual voxel edge along cell vector `axis`: |box column| / n_axis, never above the
  // nominal spacing because n_axis = ceil(|box column| / spacing).
  double h(int axis) const;
};

// Box rules (BOX keyword, see BoxMode in Options.hpp):
//   MANUAL: orthorhombic box from XLO..ZHI, origin (XLO, YLO, ZLO).
//   CELL:   the frame's cell, origin (0,0,0); error if the frame has no usable cell.
//   AUTO:   the frame's cell if usable and not more than 10x the coordinate extent of the
//           selected atoms, else a bounding box of the selected atoms with 2*spacing
//           padding, origin at their minimum corner -- byte-for-byte the pre-audit BLS path.
// n_i = max(1, ceil(|box column i| / spacing)). Fails (err set) when the two grid
// arrays would exceed 80% of available RAM.
bool deriveGridSpec(const BLSConfig& config, const Frame& frame,
                    const std::vector<int>* selection, GridSpec& spec, std::string& err);

// Configure `grid` for `spec` (always re-configures: box, origin and periodicity can
// change between frames even when the dimensions do not).
void configureGrid(Grid& grid, const GridSpec& spec);

}  // namespace bls
