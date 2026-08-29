#include "grid/Grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include "common/Types.hpp"

namespace bls {

void Grid::configure(int nx, int ny, int nz, double spacing, const Mat3& box, const Vec3& origin,
                     BoxPeriodicity periodicity) {
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  spacing_ = spacing;
  box_ = box;
  origin_ = origin;
  periodicity_ = periodicity;
  inverseBox_ = inverse(box_);
  cellVecX_ = box_.column(0) / static_cast<double>(nx_);
  cellVecY_ = box_.column(1) / static_cast<double>(ny_);
  cellVecZ_ = box_.column(2) / static_cast<double>(nz_);
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

  // Per-axis conversion from a Cartesian distance to a lattice-index reach.
  //
  // The neighbour stencil used to be sized from the single scalar
  // maxEdgeLength_ = max(|cellVecX|, |cellVecY|, |cellVecZ|), which understates
  // the number of index steps needed along any *shorter* cell vector and so
  // truncated the rasterised sphere along the short axes. Cubic grids were
  // exact (all three edges equal, and the ceil() supplied the rest), which is
  // why this went unnoticed; the reference triclinic basis lost 13.7% of the
  // voxels it should have marked at radius 5.
  //
  // The exact statement: for a voxel-cell matrix V (columns cellVecX/Y/Z), a
  // Cartesian displacement d corresponds to the index displacement V^-1 d, so
  // |(V^-1 d)_k| <= ||row_k(V^-1)|| * |d|. Since V = box * diag(1/nx,1/ny,1/nz),
  // row k of V^-1 is n_k times row k of inverseBox_ -- no second inversion, and
  // no new failure mode beyond the one inverse(box_) above already covers.
  indexReachPerLength_[0] = static_cast<double>(nx_) * norm(inverseBox_.row(0));
  indexReachPerLength_[1] = static_cast<double>(ny_) * norm(inverseBox_.row(1));
  indexReachPerLength_[2] = static_cast<double>(nz_) * norm(inverseBox_.row(2));

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
  // Fold into the cell only when the box is a periodic one. For a synthesised
  // bounding box, folding an outlying atom onto the opposite face would be a
  // fabrication; it is left where it is, and the stencil below clips whatever
  // falls outside the grid. See BoxPeriodicity in bls/Options.hpp.
  if (periodicity_ == BoxPeriodicity::Periodic) {
    frac.x -= std::floor(frac.x);
    frac.y -= std::floor(frac.y);
    frac.z -= std::floor(frac.z);
  }
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
  // One reach per axis. A single scalar reach truncates the sphere along
  // whichever cell vector is shorter than the longest; see the derivation of
  // indexReachPerLength_ in configure(). The +0.5 covers the atom's own
  // position within its base voxel, which contributes up to half an index step.
  int reach[3] = {0, 0, 0};
  if (radius > 0.0) {
    for (int k = 0; k < 3; ++k) {
      reach[k] = static_cast<int>(std::ceil(anyRadius * indexReachPerLength_[k] + 0.5));
    }
  }
  const bool periodic = periodicity_ == BoxPeriodicity::Periodic;

  auto handleAtom = [&](const Vec3& pos) {
    Vec3 wrapped;
    Vec3 frac = fractional(pos, &wrapped);
    int ix = static_cast<int>(std::floor(frac.x * nx_));
    int iy = static_cast<int>(std::floor(frac.y * ny_));
    int iz = static_cast<int>(std::floor(frac.z * nz_));
    if (periodic) {
      // frac is in [0,1) by construction, so this only catches the case where
      // a value a hair below 1 rounds up to exactly nx_.
      if (ix >= nx_) ix = nx_ - 1;
      if (iy >= ny_) iy = ny_ - 1;
      if (iz >= nz_) iz = nz_ - 1;
    }
    // Non-periodic: an atom genuinely outside the grid keeps its out-of-range
    // base index. Clamping it would teleport its footprint onto the near face;
    // visit() clips instead, so it still marks the in-range voxels that are
    // actually within reach of where it is.

    const auto visit = [&](int vx, int vy, int vz, bool isBase) {
      int sx = vx, sy = vy, sz = vz;  // storage indices
      if (periodic) {
        // Wrap the STORAGE index but measure the distance from the unwrapped
        // voxel centre, which is the correct periodic image. Wrapping first and
        // then measuring would compare against a voxel on the far side of the
        // box and reject every neighbour that crosses a face.
        sx = ((vx % nx_) + nx_) % nx_;
        sy = ((vy % ny_) + ny_) % ny_;
        sz = ((vz % nz_) + nz_) % nz_;
      } else if (vx < 0 || vy < 0 || vz < 0 || vx >= nx_ || vy >= ny_ || vz >= nz_) {
        return;
      }
      Vec3 center = voxelCenter(vx, vy, vz);
      double dist = norm(wrapped - center);
      if (mode == OccupancyMode::Any) {
        if (radius <= 0.0) {
          if (!isBase) return;
        } else if (dist > anyRadius + 1e-9) {
          return;
        }
        std::size_t idx = index(sx, sy, sz);
#ifdef _OPENMP
#pragma omp atomic write
#endif
        occ_[idx] = 1;
      } else {
        if (allRadius <= 0.0) return;
        if (dist > allRadius - 1e-9) return;
        std::size_t idx = index(sx, sy, sz);
#ifdef _OPENMP
#pragma omp atomic write
#endif
        occ_[idx] = 1;
      }
    };

    visit(ix, iy, iz, true);
    if (reach[0] > 0 || reach[1] > 0 || reach[2] > 0) {
      for (int dx = -reach[0]; dx <= reach[0]; ++dx) {
        for (int dy = -reach[1]; dy <= reach[1]; ++dy) {
          for (int dz = -reach[2]; dz <= reach[2]; ++dz) {
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
