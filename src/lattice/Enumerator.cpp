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

Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz) {
  Mat3 invBasis = inverse(basis);

  std::vector<Vec3> corners;
  corners.reserve(8);
  for (int i = 0; i <= 1; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        corners.emplace_back(static_cast<double>(i * nx), static_cast<double>(j * ny),
                             static_cast<double>(k * nz));
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

}  // namespace bls

