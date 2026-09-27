// Created: 2026-09-27T18:50+02:00 | by: CC audit (manager) | purpose: lattice-origin offset sets for bls_e7 and bls_sizefloor (audit finding S2-3, brief R9)
//
// An offset is a translation of the BLS probe lattice (deck keyword LATTICE_ORIGIN, voxels),
// not of the atoms, so the occupancy grid -- and with it the DFS ground truth -- is the same at
// every offset by construction.
//
//   halton64   the audited set. Offset 0 = (0,0,0), the historical anchoring. Offset i = 1..63
//              = B * (H2(i), H3(i), H5(i)): point i of the Halton sequence in bases 2, 3, 5 (the
//              construction E6 uses), placed in one conventional cell of the probe lattice, B = its
//              basis in voxels (Analyzer::probeBasis; the commensurate basis under PBC). 64
//              quasi-random offsets spanning the 3-D offset space (brief: ">= 64 random offsets").
//   diagonal8  the pre-audit set, kept for the PBC-none reproduction only. Offset k = 0..7 =
//              (-k,-k,-k) voxels. The pre-audit tools translated the ATOMS by (k,k,k) voxels inside
//              a fixed cell; moving the lattice by the opposite integer vector yields the same
//              partition shifted by -k voxels (llround(s - k) = llround(s) - k), so every count the
//              tools print is unchanged.
#pragma once

#include <string>

#include "common/Types.hpp"

namespace bls_tools {

inline double halton(int i, int base) {
  double f = 1.0, r = 0.0;
  while (i > 0) {
    f /= base;
    r += f * (i % base);
    i /= base;
  }
  return r;
}

// Number of offsets in `set`, or 0 for an unknown set name.
inline int offsetCount(const std::string& set) {
  if (set == "halton64") return 64;
  if (set == "diagonal8") return 8;
  return 0;
}

// LATTICE_ORIGIN (voxels) of offset `i` of `set`; `basis` = the probe basis in voxels.
inline bls::Vec3 latticeOrigin(const std::string& set, int i, const bls::Mat3& basis) {
  if (set == "diagonal8") return bls::Vec3{double(-i), double(-i), double(-i)};
  if (i == 0) return bls::Vec3{0.0, 0.0, 0.0};
  return basis * bls::Vec3{halton(i, 2), halton(i, 3), halton(i, 5)};
}

}  // namespace bls_tools
