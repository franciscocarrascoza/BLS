#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "bls/Options.hpp"
#include "common/Types.hpp"

namespace bls {

class Grid {
 public:
  void configure(int nx, int ny, int nz, double spacing, const Mat3& box, const Vec3& origin);

  void reset();

  void rasterize(const std::vector<Vec3>& positions, const std::vector<int>* selection,
                 double cutoffLength, OccupancyMode mode);

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }
  double spacing() const { return spacing_; }

  std::size_t size() const { return static_cast<std::size_t>(nx_) * ny_ * nz_; }

  const std::vector<uint8_t>& occupancy() const { return occ_; }
  std::vector<uint8_t>& occupancy() { return occ_; }

  std::vector<uint8_t>& visited() { return visited_; }
  const std::vector<uint8_t>& visited() const { return visited_; }

  const Mat3& box() const { return box_; }
  const Vec3& origin() const { return origin_; }

 private:
  // Grants tests/test_geometry.cpp access to the geometry primitives below.
  // They are pure functions of the configured box with no class invariant to
  // protect; they are private only because production code has no caller
  // outside Grid. Property-testing them directly is the whole point of the
  // accessor -- exercising them through rasterize() alone cannot separate a
  // fault in fractional() from a fault in the stencil that consumes it.
  friend struct GridTestAccess;

  std::size_t index(int ix, int iy, int iz) const;
  Vec3 fractional(const Vec3& pos, Vec3* wrapped) const;
  Vec3 voxelCenter(int ix, int iy, int iz) const;

  int nx_{0}, ny_{0}, nz_{0};
  double spacing_{1.0};
  Mat3 box_{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};
  Mat3 inverseBox_{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};
  Vec3 origin_{0.0, 0.0, 0.0};
  Vec3 cellVecX_{0.0, 0.0, 0.0};
  Vec3 cellVecY_{0.0, 0.0, 0.0};
  Vec3 cellVecZ_{0.0, 0.0, 0.0};
  double maxEdgeLength_{1.0};
  double maxCornerDistance_{0.0};

  std::vector<uint8_t> occ_;
  std::vector<uint8_t> visited_;
};

}  // namespace bls
