#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace bls {

struct SkipDFSConfig {
  int nx{0};
  int ny{0};
  int nz{0};
  int connectivity{6};
  // BLS's refinement stride (config keyword SKIP). Unrelated to
  // cluster::skipDFS's jump distance despite the shared word.
  int refinementStride{3};
  // Periodic grid (deck PBC xyz): neighbour indices wrap modulo the grid dimensions
  // instead of stopping at the faces. Default false = pre-audit behaviour.
  bool periodic{false};
};

struct SkipDFSResult {
  int componentSize{0};
  int seedHits{0};
  std::size_t refinedVoxels{0};
};

class SkipDFS {
 public:
  SkipDFS(const SkipDFSConfig& cfg, const std::vector<uint8_t>& occupancy,
          std::vector<uint8_t>& visited);

  // Walk the component containing (x,y,z). When `labels` is non-null, every
  // voxel claimed by this walk is stamped with `labelValue`; the caller owns
  // the buffer and the numbering. Null (the default) writes nothing and costs
  // one predictable branch per popped voxel.
  int runFrom(int x, int y, int z, std::vector<int>* labels = nullptr, int labelValue = -1);

  std::size_t refinedVoxels() const { return refinedVoxels_; }

 private:
  std::size_t index(int x, int y, int z) const;
  void push(int idx);
  template <bool Periodic>
  int walk(int x, int y, int z, std::vector<int>* labels, int labelValue);

  SkipDFSConfig cfg_;
  const std::vector<uint8_t>& occ_;
  std::vector<uint8_t>& visited_;
  std::vector<std::array<int, 3>> directions_;
  std::vector<int> stack_;
  std::size_t refinedVoxels_{0};
};

}  // namespace bls
