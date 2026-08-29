#include "lattice/Enumerator.hpp"

#include <algorithm>
#include <cmath>
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

void Enumerator::build(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
                       const std::vector<uint8_t>* occupancy) {
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

