#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include "common/Types.hpp"

namespace bls {

class Enumerator {
 public:
  // Enumerate all lattice points in the grid volume, independent of occupancy.
  Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz);

  // Occupancy-filtered: enumerate the lattice sites (identical set to the
  // constructor above), map each to its containing voxel, and admit it as a
  // seed iff that voxel is occupied. At most one seed per lattice site.
  Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
             const std::vector<uint8_t>& occupancy);

  struct Seed {
    int x;
    int y;
    int z;
  };

  // DEFECTIVE — retained only to reproduce pre-fix results for the record.
  // Do not use for new work.
  //
  // This was the occupancy-guided seed selection used for every BLS number
  // produced before 2026-08-29. It scanned occupied voxels and admitted each
  // one that fell within a FIXED fractional radius of 0.5 of ANY centering
  // offset. That test is not lattice membership: the per-offset spheres are
  // far larger than the Voronoi half-cell of a multi-offset lattice, so for
  // centered lattices they overlap and their union covers almost the whole
  // grid. Measured on a fully occupied 60^3 grid at dNN = 2.4 voxels, it
  // returned seed ratios of P:I:F = 1.000 : 1.783 : 1.903 against the correct
  // lattice-density ratios 1.000 : 1.299 : 1.414 — and for cubic-F, the
  // project default, 216000 seeds out of 216000 voxels. BLS as benchmarked
  // was therefore DFS with a lattice-shaped preamble.
  //
  // The lattice geometry itself (Basis.cpp computeDmin + the dNN rescaling in
  // BLS.cpp) was verified correct and is untouched; the defect was confined
  // to seed selection.
  static Enumerator legacyRadiusSelection(const Mat3& basis, const std::vector<Vec3>& offsets,
                                          int nx, int ny, int nz,
                                          const std::vector<uint8_t>& occupancy);

  template <typename F>
  void forEach(F&& func) const {
    for (const auto& seed : seeds_) {
      func(seed);
    }
  }

  int count() const { return static_cast<int>(seeds_.size()); }

 private:
  struct LegacyRadiusTag {};
  Enumerator(LegacyRadiusTag, const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
             int nz, const std::vector<uint8_t>& occupancy);

  // Shared core for both public constructors, dispatching on whether occupancy
  // is available. Without it, walks every lattice index whose site can land in
  // the grid, appends the in-range ones, then sorts and deduplicates by voxel.
  // With it, defers to buildFromOccupancy, which returns the identical set.
  void build(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
             const std::vector<uint8_t>* occupancy);

  // Occupancy-driven path taken by build() when occupancy is supplied: iterate
  // occupied voxels and ask which lattice sites round onto them, rather than
  // sweeping the whole grid volume. Produces the same seed set in the same
  // order -- see the derivation above its definition in Enumerator.cpp.
  void buildFromOccupancy(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                          int nz, const std::vector<uint8_t>& occupancy);

  std::vector<Seed> seeds_;
};

}  // namespace bls

