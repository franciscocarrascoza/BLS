#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>

#include "cluster/Algorithms.hpp"

using namespace bls;

int main() {
  // Create a simple 10x10x10 grid with some occupied voxels
  int nx = 10, ny = 10, nz = 10;
  std::size_t totalSize = nx * ny * nz;

  std::vector<uint8_t> occupancy(totalSize, 0);
  std::vector<uint8_t> visited1(totalSize, 0);
  std::vector<uint8_t> visited2(totalSize, 0);

  // Create two clusters:
  // Cluster 1: (2,2,2) to (4,4,4)
  for (int i = 2; i <= 4; ++i) {
    for (int j = 2; j <= 4; ++j) {
      for (int k = 2; k <= 4; ++k) {
        occupancy[i * ny * nz + j * nz + k] = 1;
      }
    }
  }

  // Cluster 2: (6,6,6) to (8,8,8)
  for (int i = 6; i <= 8; ++i) {
    for (int j = 6; j <= 8; ++j) {
      for (int k = 6; k <= 8; ++k) {
        occupancy[i * ny * nz + j * nz + k] = 1;
      }
    }
  }

  // Run both implementations
  auto result_basic = cc3d(nx, ny, nz, 6, occupancy, visited1);
  auto result_optimized = cc3dOptimized(nx, ny, nz, 6, occupancy, visited2);

  // Check that they produce the same results
  bool passed = true;

  if (result_basic.nclusters != result_optimized.nclusters) {
    std::cerr << "FAIL: Different cluster counts - basic: " << result_basic.nclusters
              << ", optimized: " << result_optimized.nclusters << "\n";
    passed = false;
  }

  if (result_basic.maxCluster != result_optimized.maxCluster) {
    std::cerr << "FAIL: Different max cluster sizes - basic: " << result_basic.maxCluster
              << ", optimized: " << result_optimized.maxCluster << "\n";
    passed = false;
  }

  if (result_basic.visitedVoxels != result_optimized.visitedVoxels) {
    std::cerr << "FAIL: Different visited voxel counts - basic: " << result_basic.visitedVoxels
              << ", optimized: " << result_optimized.visitedVoxels << "\n";
    passed = false;
  }

  // Check cluster sizes (should be same after sorting)
  auto sizes_basic = result_basic.clusterSizes;
  auto sizes_optimized = result_optimized.clusterSizes;
  std::sort(sizes_basic.begin(), sizes_basic.end());
  std::sort(sizes_optimized.begin(), sizes_optimized.end());

  if (sizes_basic != sizes_optimized) {
    std::cerr << "FAIL: Different cluster size distributions\n";
    passed = false;
  }

  // Expected results for this test case
  if (result_basic.nclusters != 2) {
    std::cerr << "FAIL: Expected 2 clusters, got " << result_basic.nclusters << "\n";
    passed = false;
  }

  if (result_basic.maxCluster != 27) {  // 3x3x3 = 27 voxels per cluster
    std::cerr << "FAIL: Expected max cluster size 27, got " << result_basic.maxCluster << "\n";
    passed = false;
  }

  // Performance note: basic should be slower but produce same results
  if (result_basic.elapsedMs == 0.0 || result_optimized.elapsedMs == 0.0) {
    std::cerr << "WARNING: One implementation took <0.001ms, may need larger test\n";
  }

  if (passed) {
    std::cout << "✓ CC3D fairness test PASSED\n";
    std::cout << "  Both implementations found " << result_basic.nclusters << " clusters\n";
    std::cout << "  Max cluster size: " << result_basic.maxCluster << " voxels\n";
    std::cout << "  Timing - basic: " << result_basic.elapsedMs << " ms, "
              << "optimized: " << result_optimized.elapsedMs << " ms\n";
    std::cout << "  Correctness verified: identical cluster results\n";
    return 0;
  } else {
    std::cerr << "✗ CC3D fairness test FAILED\n";
    return 1;
  }
}
