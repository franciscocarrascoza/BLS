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

  // BLS seeding as the manuscript describes it (§2, §2.1.2; audit decision D3): every
  // lattice probe site whose rounded voxel lies in the grid is evaluated once -- round
  // to the voxel, read its occupancy, O(1) -- and the occupied ones are the seed set S,
  // put in lexicographic (x,y,z) order and de-duplicated (two sites can round onto one
  // voxel when dNN < sqrt(3) voxels). probes() is the number of sites evaluated, the m
  // of §2.1.2. The seed set and its order are those of occupiedSweep() below; only the
  // work differs (m site evaluations instead of a sweep over all M^3 voxels).
  // `origin` (voxels, deck LATTICE_ORIGIN) translates every site before rounding:
  // site = basis * (idx + offset) + origin. The default zero leaves the arithmetic exact.
  Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
             const std::vector<uint8_t>& occupancy, const Vec3& origin = Vec3{});

  // Periodic grid (deck PBC xyz; audit decision D10). `basis` must be a commensurate
  // basis from commensurateCubicBasis() (Basis.hpp): diagonal, a_i = N_i / n_i with
  // n_i integral, so the lattice tiles the periodic grid. Sites idx + offset with
  // idx_i in [0, n_i) are evaluated; voxel = llround(site) wrapped modulo N_i.
  // probes() = n_x * n_y * n_z * offsets.size(). Throws std::invalid_argument for a
  // basis that is not commensurate with the grid.
  struct PeriodicTag {};
  Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
             const std::vector<uint8_t>& occupancy, PeriodicTag, const Vec3& origin = Vec3{});

  // BLS seeding before audit decision D3, kept as a test oracle only: sweeps every voxel
  // (word-skipping empty runs) and tests each occupied one for lattice membership. Same
  // seed set in the same order as the probe evaluation above, by construction (see the
  // derivation in Enumerator.cpp); probes() is 0 because no site is evaluated.
  static Enumerator occupiedSweep(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                                  int ny, int nz, const std::vector<uint8_t>& occupancy,
                                  const Vec3& origin = Vec3{});

  // Test oracles for the separable probe loop (audit 28.09.26, BLS-PERF2): the same probe
  // evaluation with that shortcut disabled -- the row-by-row loops that ran before it, which
  // round every site's coordinates. Same seeds, same order and same probes() as the
  // constructors above, by the argument in build().
  static Enumerator rowWise(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                            int nz, const std::vector<uint8_t>& occupancy,
                            const Vec3& origin = Vec3{});
  static Enumerator rowWisePeriodic(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                                    int ny, int nz, const std::vector<uint8_t>& occupancy,
                                    const Vec3& origin = Vec3{});

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

  // Lattice probe sites evaluated while building (m of §2.1.2); 0 for occupiedSweep().
  long long probes() const { return probes_; }

 private:
  Enumerator() = default;
  struct LegacyRadiusTag {};
  Enumerator(LegacyRadiusTag, const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
             int nz, const std::vector<uint8_t>& occupancy);

  // Probe evaluation shared by the non-periodic constructors: walks every lattice
  // index whose site can land in the grid, evaluates the in-range sites (counting
  // them in probes_), keeps those on an occupied voxel (all of them when occupancy is
  // null), then sorts and de-duplicates by voxel.
  // A diagonal basis (every cubic lattice) takes the separable path: each coordinate rounded
  // once per index along its own axis, the sites the product of the three lists, the
  // per-offset runs merged rather than sorted. allowSeparable = false is rowWise().
  void build(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
             const std::vector<uint8_t>* occupancy, const Vec3& origin = Vec3{},
             bool allowSeparable = true);

  // Periodic probe evaluation (PeriodicTag constructor).
  void buildPeriodic(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
                     const std::vector<uint8_t>& occupancy, const Vec3& origin,
                     bool allowSeparable = true);

  // occupiedSweep(): iterate occupied voxels and ask which lattice sites round onto
  // them. Produces the same seed set in the same order as build() -- see the
  // derivation above its definition in Enumerator.cpp.
  void buildFromOccupancy(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                          int nz, const std::vector<uint8_t>& occupancy, const Vec3& origin);

  std::vector<Seed> seeds_;
  long long probes_{0};
};

}  // namespace bls

