#include "lattice/Enumerator.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "common/Types.hpp"

namespace bls {

namespace {

struct Bounds {
  double minX = std::numeric_limits<double>::infinity();
  double maxX = -std::numeric_limits<double>::infinity();
  double minY = std::numeric_limits<double>::infinity();
  double maxY = -std::numeric_limits<double>::infinity();
  double minZ = std::numeric_limits<double>::infinity();
  double maxZ = -std::numeric_limits<double>::infinity();
};

}  // namespace

// Occupancy-driven enumeration. Same seed set as the volume sweep below, found
// by inverting the question instead of sweeping: for each OCCUPIED voxel, is
// there a lattice site that rounds to it?
//
// Why it is worth having: the volume sweep visits every lattice site in the
// grid -- ~2.3M on E1/Ic -- to find the ~2.2k that land on the ~21k occupied
// voxels. E1 occupancy is 0.09%, so >99.9% of that work is discarded, and
// enumeration was 69-77% of BLS's total pipeline time on all four E1 systems.
//
// Why the seed set is provably identical, not merely equal in testing:
// llround(site) == v implies |site_k - v_k| <= 0.5 for every k. Writing
// idx = invBasis*site - offset and f0 = invBasis*v - offset gives
// idx - f0 = invBasis*(site - v), so
//     |idx_k - f0_k| <= 0.5 * sum_j |invBasis[k][j]|
// which is the per-axis search half-width below. Every lattice index that
// could round into v therefore lies inside the searched box; nothing outside
// it needs testing. The acceptance test itself is byte-for-byte the volume
// sweep's -- site = basis * (latticeIdx + offset), then llround -- so a
// candidate is admitted on exactly the same arithmetic, with no tolerance and
// no re-derived predicate. The bound is a superset, so widening it (the eps
// below) can only add candidates that the exact test then rejects; it can
// never change the answer.
//
// Emission order needs no sort: voxels are visited in x-major order, which is
// the lexicographic (x,y,z) order the volume path sorts into, and each voxel is
// tested once, which is what its dedup achieves.
void Enumerator::buildFromOccupancy(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                                    int ny, int nz, const std::vector<uint8_t>& occupancy) {
  const Mat3 invBasis = inverse(basis);

  // Half-width of the candidate box, per lattice-index axis. Derived above; the
  // epsilon covers rounding in the bound's own arithmetic only.
  const double eps = 1e-9;
  double half[3];
  for (int k = 0; k < 3; ++k) {
    Vec3 r = invBasis.row(k);
    half[k] = 0.5 * (std::fabs(r.x) + std::fabs(r.y) + std::fabs(r.z)) + eps;
  }

  std::vector<Vec3> offsetsReal;
  offsetsReal.reserve(offsets.size());
  for (const auto& o : offsets) offsetsReal.push_back(basis * o);

  const std::size_t total = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                            static_cast<std::size_t>(nz);
  const std::size_t planeYZ = static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);

  for (std::size_t base = 0; base < total; ) {
    // Occupancy is a few tenths of a percent dense here, so skip runs of empty
    // voxels a machine word at a time rather than a byte at a time.
    if (base + 8 <= total && (base & 7u) == 0) {
      std::uint64_t word;
      std::memcpy(&word, occupancy.data() + base, 8);
      if (word == 0) { base += 8; continue; }
    }
    if (!occupancy[base]) { ++base; continue; }

    const int vx = static_cast<int>(base / planeYZ);
    const std::size_t rem = base - static_cast<std::size_t>(vx) * planeYZ;
    const int vy = static_cast<int>(rem / static_cast<std::size_t>(nz));
    const int vz = static_cast<int>(rem % static_cast<std::size_t>(nz));
    const Vec3 v{static_cast<double>(vx), static_cast<double>(vy), static_cast<double>(vz)};

    bool admitted = false;
    for (std::size_t oi = 0; oi < offsets.size() && !admitted; ++oi) {
      const Vec3 f0 = invBasis * (v - offsetsReal[oi]);
      const int ixLo = static_cast<int>(std::ceil(f0.x - half[0]));
      const int ixHi = static_cast<int>(std::floor(f0.x + half[0]));
      if (ixLo > ixHi) continue;
      const int iyLo = static_cast<int>(std::ceil(f0.y - half[1]));
      const int iyHi = static_cast<int>(std::floor(f0.y + half[1]));
      if (iyLo > iyHi) continue;
      const int izLo = static_cast<int>(std::ceil(f0.z - half[2]));
      const int izHi = static_cast<int>(std::floor(f0.z + half[2]));
      if (izLo > izHi) continue;

      for (int ix = ixLo; ix <= ixHi && !admitted; ++ix) {
        for (int iy = iyLo; iy <= iyHi && !admitted; ++iy) {
          for (int iz = izLo; iz <= izHi; ++iz) {
            Vec3 latticeIdx{static_cast<double>(ix), static_cast<double>(iy),
                            static_cast<double>(iz)};
            Vec3 site = basis * (latticeIdx + offsets[oi]);
            if (static_cast<int>(std::llround(site.x)) == vx &&
                static_cast<int>(std::llround(site.y)) == vy &&
                static_cast<int>(std::llround(site.z)) == vz) {
              admitted = true;
              break;
            }
          }
        }
      }
    }
    if (admitted) seeds_.push_back(Seed{vx, vy, vz});
    ++base;
  }
}

void Enumerator::build(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
                       const std::vector<uint8_t>* occupancy) {
  if (occupancy) {
    buildFromOccupancy(basis, offsets, nx, ny, nz, *occupancy);
    return;
  }

  Mat3 invBasis = inverse(basis);

  // A site s is kept iff llround(s_k) is in [0, n_k), i.e. iff s lies in the
  // half-open box [-0.5, n_k - 0.5). Bound the lattice indices with the corners
  // of THAT box, not of [0, n_k] -- otherwise the low-side half-voxel shell is
  // outside the hull and has to be recovered by padding. With the corners
  // stated exactly, a padding of 1 covers only the floor/ceil of the hull
  // itself and is not load-bearing: widening it further adds no sites.
  std::vector<Vec3> corners;
  corners.reserve(8);
  for (int i = 0; i <= 1; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        corners.emplace_back(i ? static_cast<double>(nx) - 0.5 : -0.5,
                             j ? static_cast<double>(ny) - 0.5 : -0.5,
                             k ? static_cast<double>(nz) - 0.5 : -0.5);
      }
    }
  }

  std::vector<Seed> provisional;

  for (const auto& offset : offsets) {
    Vec3 offsetReal = basis * offset;
    Bounds bnd;
    for (const auto& corner : corners) {
      Vec3 idx = invBasis * (corner - offsetReal);
      bnd.minX = std::min(bnd.minX, idx.x);
      bnd.maxX = std::max(bnd.maxX, idx.x);
      bnd.minY = std::min(bnd.minY, idx.y);
      bnd.maxY = std::max(bnd.maxY, idx.y);
      bnd.minZ = std::min(bnd.minZ, idx.z);
      bnd.maxZ = std::max(bnd.maxZ, idx.z);
    }

    int ixMin = static_cast<int>(std::floor(bnd.minX)) - 1;
    int ixMax = static_cast<int>(std::ceil(bnd.maxX)) + 1;
    int iyMin = static_cast<int>(std::floor(bnd.minY)) - 1;
    int iyMax = static_cast<int>(std::ceil(bnd.maxY)) + 1;
    int izMin = static_cast<int>(std::floor(bnd.minZ)) - 1;
    int izMax = static_cast<int>(std::ceil(bnd.maxZ)) + 1;

    for (int ix = ixMin; ix <= ixMax; ++ix) {
      for (int iy = iyMin; iy <= iyMax; ++iy) {
        for (int iz = izMin; iz <= izMax; ++iz) {
          Vec3 latticeIdx{static_cast<double>(ix), static_cast<double>(iy),
                          static_cast<double>(iz)};
          Vec3 site = basis * (latticeIdx + offset);
          int vx = static_cast<int>(std::llround(site.x));
          int vy = static_cast<int>(std::llround(site.y));
          int vz = static_cast<int>(std::llround(site.z));
          if (vx < 0 || vy < 0 || vz < 0 || vx >= nx || vy >= ny || vz >= nz) continue;
          if (occupancy) {
            // Reject here rather than after the sort: on sparse grids this is
            // the difference between materialising every lattice site in the
            // box and materialising only the seeds.
            std::size_t idx = static_cast<std::size_t>(vx) * static_cast<std::size_t>(ny) *
                                  static_cast<std::size_t>(nz) +
                              static_cast<std::size_t>(vy) * static_cast<std::size_t>(nz) +
                              static_cast<std::size_t>(vz);
            if (!(*occupancy)[idx]) continue;
          }
          provisional.push_back(Seed{vx, vy, vz});
        }
      }
    }
  }

  std::sort(provisional.begin(), provisional.end(),
            [](const Seed& a, const Seed& b) {
              if (a.x != b.x) return a.x < b.x;
              if (a.y != b.y) return a.y < b.y;
              return a.z < b.z;
            });
  provisional.erase(std::unique(provisional.begin(), provisional.end(),
                                [](const Seed& a, const Seed& b) {
                                  return a.x == b.x && a.y == b.y && a.z == b.z;
                                }),
                    provisional.end());

  seeds_ = std::move(provisional);
}

Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz) {
  build(basis, offsets, nx, ny, nz, nullptr);
}

// Same lattice sites as the constructor above, restricted to those whose voxel
// is occupied. Seed selection is decided entirely by lattice geometry;
// occupancy only filters.
Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz, const std::vector<uint8_t>& occupancy) {
  build(basis, offsets, nx, ny, nz, &occupancy);
}

Enumerator Enumerator::legacyRadiusSelection(const Mat3& basis, const std::vector<Vec3>& offsets,
                                             int nx, int ny, int nz,
                                             const std::vector<uint8_t>& occupancy) {
  return Enumerator(LegacyRadiusTag{}, basis, offsets, nx, ny, nz, occupancy);
}

// DEFECTIVE — see the note on Enumerator::legacyRadiusSelection in Enumerator.hpp.
// Preserved verbatim so pre-2026-08-29 results can be reproduced. The bug is the
// unconditional per-offset radius test below: `< 0.25` is a fixed fractional
// radius of 0.5 applied once per centering offset, so multi-offset centerings
// capture overlapping spheres covering nearly the entire grid rather than the
// lattice sites. Not to be repaired by retuning the constant — a radius fitted
// to cubic P/I/F would still be wrong for hexagonal and triclinic.
Enumerator::Enumerator(LegacyRadiusTag, const Mat3& basis, const std::vector<Vec3>& offsets,
                       int nx, int ny, int nz, const std::vector<uint8_t>& occupancy) {
  // Frozen copy of the pre-2026-08-29 bls::inverse(), which returned the
  // TRANSPOSED inverse. This constructor exists solely to reproduce pre-fix
  // numbers, so it must reproduce the pre-fix arithmetic exactly -- not merely
  // the pre-fix algebra. Calling the corrected inverse() and transposing it
  // back is not good enough: the radius test below is a strict `< 0.25` on a
  // knife edge that voxels land on exactly (at dNN = 2.4 voxels, every voxel
  // 1.2 from a site), so a one-ulp difference in the matrix flips hundreds of
  // seeds. Do not replace this with a call to inverse().
  auto frozenPreFixInverse = [](const Mat3& m) {
    Vec3 c0 = cross(m.cols[1], m.cols[2]);
    Vec3 c1 = cross(m.cols[2], m.cols[0]);
    Vec3 c2 = cross(m.cols[0], m.cols[1]);
    double det = dot(m.cols[0], c0);
    return Mat3{c0 / det, c1 / det, c2 / det};
  };
  Mat3 invBasis = frozenPreFixInverse(basis);

  // Precompute real-space offsets
  std::vector<Vec3> offsetsReal;
  offsetsReal.reserve(offsets.size());
  for (const auto& o : offsets) {
    offsetsReal.push_back(basis * o);
  }

  // Scan the occupancy grid sequentially (cache-friendly).
  // For each occupied voxel, check if it lies on any lattice site.
  for (int x = 0; x < nx; ++x) {
    std::size_t xStride = static_cast<std::size_t>(x) * static_cast<std::size_t>(ny) *
                          static_cast<std::size_t>(nz);
    for (int y = 0; y < ny; ++y) {
      std::size_t xyStride = xStride + static_cast<std::size_t>(y) * static_cast<std::size_t>(nz);
      for (int z = 0; z < nz; ++z) {
        if (!occupancy[xyStride + static_cast<std::size_t>(z)]) continue;

        Vec3 pos{static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)};
        for (const auto& oReal : offsetsReal) {
          Vec3 frac = invBasis * (pos - oReal);
          double dx = frac.x - std::round(frac.x);
          double dy = frac.y - std::round(frac.y);
          double dz = frac.z - std::round(frac.z);
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            seeds_.push_back(Seed{x, y, z});
            break;  // one match is enough
          }
        }
      }
    }
  }
}

}  // namespace bls

