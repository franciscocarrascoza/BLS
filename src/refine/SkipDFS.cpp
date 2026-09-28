#include "refine/SkipDFS.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace bls {

namespace {

std::vector<std::array<int, 3>> buildDirections(int connectivity) {
  std::vector<std::array<int, 3>> dirs;
  if (connectivity == 6) {
    dirs = {{{1, 0, 0}},
            {{-1, 0, 0}},
            {{0, 1, 0}},
            {{0, -1, 0}},
            {{0, 0, 1}},
            {{0, 0, -1}}};
  } else if (connectivity == 18) {
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          int manhattan = std::abs(dx) + std::abs(dy) + std::abs(dz);
          if (manhattan >= 1 && manhattan <= 2) {
            dirs.push_back({dx, dy, dz});
          }
        }
      }
    }
  } else {
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          if (dx == 0 && dy == 0 && dz == 0) continue;
          dirs.push_back({dx, dy, dz});
        }
      }
    }
  }
  return dirs;
}

// Periodic index wrap, identical to ((a % n) + n) % n for every int a, without the two divisions
// when a is already inside [0, n) -- the common case. The same helper as the comparison methods'
// (cluster/Algorithms.cpp, finding PBC-PERF), which BLS's own refinement did not get then
// (audit 28.09.26, BLS-PERF2).
inline int wrapIndex(int a, int n) {
  if (static_cast<unsigned>(a) < static_cast<unsigned>(n)) return a;
  const int r = a % n;
  return r < 0 ? r + n : r;
}

}  // namespace

SkipDFS::SkipDFS(const SkipDFSConfig& cfg, const std::vector<uint8_t>& occupancy,
                 std::vector<uint8_t>& visited)
    : cfg_(cfg), occ_(occupancy), visited_(visited), directions_(buildDirections(cfg.connectivity)) {
  stack_.reserve(1024);
}

std::size_t SkipDFS::index(int x, int y, int z) const {
  return static_cast<std::size_t>(x) * cfg_.ny * cfg_.nz + static_cast<std::size_t>(y) * cfg_.nz +
         static_cast<std::size_t>(z);
}

int SkipDFS::runFrom(int x, int y, int z, std::vector<int>* labels, int labelValue) {
  // The periodic walk is a separate instantiation, so PBC none runs exactly the
  // pre-audit loop (same code, same cost) and only PBC xyz pays for the wrap.
  return cfg_.periodic ? walk<true>(x, y, z, labels, labelValue)
                       : walk<false>(x, y, z, labels, labelValue);
}

template <bool Periodic>
int SkipDFS::walk(int x, int y, int z, std::vector<int>* labels, int labelValue) {
  if (x < 0 || y < 0 || z < 0 || x >= cfg_.nx || y >= cfg_.ny || z >= cfg_.nz) return 0;
  std::size_t start = index(x, y, z);
  if (!occ_[start] || visited_[start]) return 0;

  int clusterSize = 0;
  int maxSkip = std::max(1, cfg_.refinementStride);
  refinedVoxels_ = 0;

  stack_.clear();
  stack_.push_back({x, y, z});
  visited_[start] = 1;

  while (!stack_.empty()) {
    const Cell c = stack_.back();
    stack_.pop_back();
    if (labels) (*labels)[index(c.x, c.y, c.z)] = labelValue;
    ++clusterSize;
    ++refinedVoxels_;

    for (const auto& dir : directions_) {
      for (int step = 1; step <= maxSkip; ++step) {
        int nx = c.x + dir[0] * step;
        int ny = c.y + dir[1] * step;
        int nz = c.z + dir[2] * step;
        if constexpr (Periodic) {
          // |dir * step| can exceed a small grid dimension, so the wrap is a true modulo
          // outside [0, n).
          nx = wrapIndex(nx, cfg_.nx);
          ny = wrapIndex(ny, cfg_.ny);
          nz = wrapIndex(nz, cfg_.nz);
        } else {
          if (nx < 0 || ny < 0 || nz < 0 || nx >= cfg_.nx || ny >= cfg_.ny || nz >= cfg_.nz) break;
        }
        const std::size_t nidx = index(nx, ny, nz);
        if (!occ_[nidx]) break;
        if (!visited_[nidx]) {
          visited_[nidx] = 1;
          stack_.push_back({nx, ny, nz});
        }
      }
    }
  }

  return clusterSize;
}

}  // namespace bls
