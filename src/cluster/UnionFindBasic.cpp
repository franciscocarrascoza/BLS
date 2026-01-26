#include "cluster/UnionFindBasic.hpp"

namespace bls {

UnionFindBasic::UnionFindBasic(std::size_t size) : parent_(size) {
  reset();
}

void UnionFindBasic::reset() {
  for (std::size_t i = 0; i < parent_.size(); ++i) {
    parent_[i] = static_cast<int>(i);
  }
}

int UnionFindBasic::find(int x) {
  // Basic find WITHOUT path compression
  // Follow parent pointers until we reach root
  while (parent_[x] != x) {
    x = parent_[x];
  }
  return x;
}

void UnionFindBasic::unite(int x, int y) {
  // Basic union WITHOUT union-by-rank
  // Always make root of y point to root of x
  int rootX = find(x);
  int rootY = find(y);
  if (rootX != rootY) {
    parent_[rootY] = rootX;
  }
}

}  // namespace bls
