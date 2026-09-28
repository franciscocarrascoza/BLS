#include "refine/StandardDFS.hpp"

#include <cstdlib>

namespace bls {

namespace {

std::vector<std::array<int, 3>> buildDirections(int connectivity) {
  if (connectivity == 6) {
    // The traditional_dfs baseline's deltas6, in its order (cluster/Algorithms.cpp).
    return {{{-1, 0, 0}}, {{1, 0, 0}}, {{0, -1, 0}}, {{0, 1, 0}}, {{0, 0, -1}}, {{0, 0, 1}}};
  }
  std::vector<std::array<int, 3>> dirs;
  const int maxManhattan = connectivity == 18 ? 2 : 3;
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dz = -1; dz <= 1; ++dz) {
        const int manhattan = std::abs(dx) + std::abs(dy) + std::abs(dz);
        if (manhattan >= 1 && manhattan <= maxManhattan) dirs.push_back({dx, dy, dz});
      }
    }
  }
  return dirs;
}

}  // namespace

StandardDFS::StandardDFS(const SkipDFSConfig& cfg, const std::vector<uint8_t>& occupancy,
                         std::vector<uint8_t>& visited)
    : cfg_(cfg), occ_(occupancy), visited_(visited), directions_(buildDirections(cfg.connectivity)) {
  stack_.reserve(10000);  // as the baseline
}

int StandardDFS::runFrom(int x, int y, int z, std::vector<int>* labels, int labelValue) {
  return cfg_.periodic ? walk<true>(x, y, z, labels, labelValue)
                       : walk<false>(x, y, z, labels, labelValue);
}

template <bool Periodic>
int StandardDFS::walk(int x, int y, int z, std::vector<int>* labels, int labelValue) {
  const int nx = cfg_.nx, ny = cfg_.ny, nz = cfg_.nz;
  refinedVoxels_ = 0;
  if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) return 0;
  const std::size_t start = (static_cast<std::size_t>(x) * ny + y) * nz + z;
  if (!occ_[start] || visited_[start]) return 0;

  int clusterSize = 0;
  stack_.clear();
  stack_.push_back({x, y, z});
  while (!stack_.empty()) {
    const Cell c = stack_.back();
    stack_.pop_back();
    if (c.i < 0 || c.i >= nx || c.j < 0 || c.j >= ny || c.k < 0 || c.k >= nz) continue;
    const std::size_t index = (static_cast<std::size_t>(c.i) * ny + c.j) * nz + c.k;
    if (occ_[index] == 0 || visited_[index]) continue;
    visited_[index] = 1;
    if (labels) (*labels)[index] = labelValue;
    ++clusterSize;
    for (const auto& d : directions_) {
      int ni = c.i + d[0], nj = c.j + d[1], nk = c.k + d[2];
      if constexpr (Periodic) {
        ni = ni < 0 ? ni + nx : (ni >= nx ? ni - nx : ni);
        nj = nj < 0 ? nj + ny : (nj >= ny ? nj - ny : nj);
        nk = nk < 0 ? nk + nz : (nk >= nz ? nk - nz : nk);
        stack_.push_back({ni, nj, nk});
      } else {
        if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
          stack_.push_back({ni, nj, nk});
        }
      }
    }
  }
  refinedVoxels_ = static_cast<std::size_t>(clusterSize);
  return clusterSize;
}

}  // namespace bls
