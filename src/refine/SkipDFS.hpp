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
  int skip{3};
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

  int runFrom(int x, int y, int z);

  std::size_t refinedVoxels() const { return refinedVoxels_; }

 private:
  std::size_t index(int x, int y, int z) const;
  void push(int idx);

  SkipDFSConfig cfg_;
  const std::vector<uint8_t>& occ_;
  std::vector<uint8_t>& visited_;
  std::vector<std::array<int, 3>> directions_;
  std::vector<int> stack_;
  std::size_t refinedVoxels_{0};
};

}  // namespace bls
