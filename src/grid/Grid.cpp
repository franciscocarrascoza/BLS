#include "grid/Grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include "common/Types.hpp"

namespace bls {

void Grid::configure(int nx, int ny, int nz, double spacing, const Mat3& box, const Vec3& origin) {
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  spacing_ = spacing;
  box_ = box;
  origin_ = origin;
  inverseBox_ = inverse(box_);
  cellVecX_ = box_.column(0) / static_cast<double>(nx_);
  cellVecY_ = box_.column(1) / static_cast<double>(ny_);
  cellVecZ_ = box_.column(2) / static_cast<double>(nz_);
  maxEdgeLength_ =
      std::max({norm(cellVecX_), norm(cellVecY_), norm(cellVecZ_), std::numeric_limits<double>::min()});

  Vec3 halfX = cellVecX_ * 0.5;
  Vec3 halfY = cellVecY_ * 0.5;
  Vec3 halfZ = cellVecZ_ * 0.5;
  maxCornerDistance_ = 0.0;
  for (int sx : {-1, 1}) {
    for (int sy : {-1, 1}) {
      for (int sz : {-1, 1}) {
        Vec3 corner = halfX * static_cast<double>(sx) + halfY * static_cast<double>(sy) +
                      halfZ * static_cast<double>(sz);
        maxCornerDistance_ = std::max(maxCornerDistance_, norm(corner));
      }
    }
  }

  occ_.assign(size(), 0);
  visited_.assign(size(), 0);
}

void Grid::reset() {
  std::fill(occ_.begin(), occ_.end(), 0);
  std::fill(visited_.begin(), visited_.end(), 0);
}

std::size_t Grid::index(int ix, int iy, int iz) const {
  return static_cast<std::size_t>(ix) * static_cast<std::size_t>(ny_) * nz_ +
         static_cast<std::size_t>(iy) * static_cast<std::size_t>(nz_) +
         static_cast<std::size_t>(iz);
}

Vec3 Grid::fractional(const Vec3& pos, Vec3* wrapped) const {
  Vec3 rel = pos - origin_;
  Vec3 frac = inverseBox_ * rel;
  frac.x -= std::floor(frac.x);
  frac.y -= std::floor(frac.y);
  frac.z -= std::floor(frac.z);
  if (wrapped) {
    *wrapped = origin_ + box_ * frac;
  }
  return frac;
}

Vec3 Grid::voxelCenter(int ix, int iy, int iz) const {
  Vec3 frac{(static_cast<double>(ix) + 0.5) / static_cast<double>(nx_),
            (static_cast<double>(iy) + 0.5) / static_cast<double>(ny_),
            (static_cast<double>(iz) + 0.5) / static_cast<double>(nz_)};
  return origin_ + box_ * frac;
}

void Grid::rasterize(const std::vector<Vec3>& positions, const std::vector<int>* selection,
                     double cutoffLength, OccupancyMode mode) {
  reset();
  const int selectionCount =
      selection ? static_cast<int>(selection->size()) : static_cast<int>(positions.size());
  if (selectionCount == 0) return;

  const double radius = std::max(0.0, cutoffLength);
  const double anyRadius = radius + maxCornerDistance_;
  const double allRadius = radius - maxCornerDistance_;
  const int reach =
      radius > 0.0 ? static_cast<int>(std::ceil((radius + maxCornerDistance_) / maxEdgeLength_)) : 0;

  auto handleAtom = [&](const Vec3& pos) {
    Vec3 wrapped;
    Vec3 frac = fractional(pos, &wrapped);
    int ix = static_cast<int>(std::floor(frac.x * nx_));
    int iy = static_cast<int>(std::floor(frac.y * ny_));
    int iz = static_cast<int>(std::floor(frac.z * nz_));
    if (ix >= nx_) ix = nx_ - 1;
    if (iy >= ny_) iy = ny_ - 1;
    if (iz >= nz_) iz = nz_ - 1;

    const auto visit = [&](int vx, int vy, int vz, bool isBase) {
      if (vx < 0 || vy < 0 || vz < 0 || vx >= nx_ || vy >= ny_ || vz >= nz_) return;
      Vec3 center = voxelCenter(vx, vy, vz);
      double dist = norm(wrapped - center);
      if (mode == OccupancyMode::Any) {
        if (radius <= 0.0) {
          if (!isBase) return;
        } else if (dist > anyRadius + 1e-9) {
          return;
        }
        std::size_t idx = index(vx, vy, vz);
#ifdef _OPENMP
#pragma omp atomic write
#endif
        occ_[idx] = 1;
      } else {
        if (allRadius <= 0.0) return;
        if (dist > allRadius - 1e-9) return;
        std::size_t idx = index(vx, vy, vz);
#ifdef _OPENMP
#pragma omp atomic write
#endif
        occ_[idx] = 1;
      }
    };

    visit(ix, iy, iz, true);
    if (reach > 0) {
      for (int dx = -reach; dx <= reach; ++dx) {
        for (int dy = -reach; dy <= reach; ++dy) {
          for (int dz = -reach; dz <= reach; ++dz) {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            visit(ix + dx, iy + dy, iz + dz, false);
          }
        }
      }
    }
  };

  if (selection) {
    const auto& sel = *selection;
    std::size_t nsel = sel.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nsel > 256)
#endif
    for (std::size_t i = 0; i < nsel; ++i) {
      int idx = sel[i];
      if (idx < 0 || idx >= static_cast<int>(positions.size())) continue;
      handleAtom(positions[static_cast<std::size_t>(idx)]);
    }
  } else {
    std::size_t npos = positions.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(npos > 256)
#endif
    for (std::size_t i = 0; i < npos; ++i) {
      handleAtom(positions[i]);
    }
  }
}

}  // namespace bls
