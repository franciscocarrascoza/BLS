#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "refine/SkipDFS.hpp"

namespace bls {

// BLS refinement by the traditional_dfs baseline's own flood fill (deck REFINEMENT DFS; audit
// 28.09.26). The walk is the baseline's (cluster/Algorithms.cpp, traditionalDFSImpl): a stack
// of voxel coordinates, each popped voxel checked for range, occupancy and visited, then
// marked, and every neighbour pushed (in range, or wrapped under PBC). The only difference
// from the baseline is where the walks start: at BLS's seeds, not at every voxel of a raster
// scan. At CONNECTIVITY 6 the neighbour order is the baseline's (deltas6); 18 and 26 use the
// same direction sets as SkipDFS. Same interface as SkipDFS; the stride is ignored.
class StandardDFS {
 public:
  StandardDFS(const SkipDFSConfig& cfg, const std::vector<uint8_t>& occupancy,
              std::vector<uint8_t>& visited);

  // Walk the component containing (x,y,z); 0 when the voxel is out of range, empty or
  // already visited. Labels as SkipDFS::runFrom.
  int runFrom(int x, int y, int z, std::vector<int>* labels = nullptr, int labelValue = -1);

  std::size_t refinedVoxels() const { return refinedVoxels_; }

 private:
  struct Cell {
    int i, j, k;
  };
  template <bool Periodic>
  int walk(int x, int y, int z, std::vector<int>* labels, int labelValue);

  SkipDFSConfig cfg_;
  const std::vector<uint8_t>& occ_;
  std::vector<uint8_t>& visited_;
  std::vector<std::array<int, 3>> directions_;
  std::vector<Cell> stack_;
  std::size_t refinedVoxels_{0};
};

}  // namespace bls
