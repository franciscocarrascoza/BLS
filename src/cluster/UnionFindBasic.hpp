#pragma once

#include <cstddef>
#include <vector>

namespace bls {

// Basic Union-Find structure WITHOUT path compression or union-by-rank
// This implementation matches BLS's optimization level for fair benchmarking
class UnionFindBasic {
 public:
  explicit UnionFindBasic(std::size_t size);

  // Find root WITHOUT path compression
  // Time complexity: O(tree height), worst case O(n)
  int find(int x);

  // Union WITHOUT union-by-rank
  // Always makes the second root a child of the first
  void unite(int x, int y);

  // Reset the structure
  void reset();

  std::size_t size() const { return parent_.size(); }

 private:
  std::vector<int> parent_;
};

}  // namespace bls
