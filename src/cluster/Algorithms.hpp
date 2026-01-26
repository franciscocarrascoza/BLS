#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace bls {

// Enumeration of available clustering algorithms
enum class ClusterAlgorithm {
  BLS,           // Bravais Lattice Sampling (default)
  TraditionalDFS,// Traditional DFS
  SkipDFS,       // Skip-DFS without lattice sampling
  DBSCAN,        // Density-based clustering
  Hierarchical,  // Single-linkage hierarchical
  KMeans,        // K-means clustering
  Spectral,      // Simplified spectral clustering
  GCBD,          // Union-Find based clustering
  HDBSCAN,       // Hierarchical DBSCAN
  CC3D,          // Connected Components 3D (fair: basic Union-Find)
  CC3DOptimized  // Connected Components 3D (optimized: path compression + union-by-rank)
};

// Convert string to algorithm enum
ClusterAlgorithm parseAlgorithm(const std::string& name);

// Convert algorithm enum to string
std::string algorithmToString(ClusterAlgorithm algo);

// List all available algorithms
std::vector<std::string> listAlgorithms();

// Result structure for clustering operations
struct ClusterResult {
  int nclusters{0};
  int maxCluster{0};
  std::size_t visitedVoxels{0};
  double elapsedMs{0.0};
  std::vector<int> clusterSizes;
};

// Parameters for clustering algorithms
struct ClusterParams {
  int nx{0}, ny{0}, nz{0};     // Grid dimensions
  int skip{3};                  // Skip distance for Skip-DFS
  double eps{3.0};              // Epsilon for DBSCAN
  int minPts{10};               // MinPts for DBSCAN
  int k{20};                    // Number of clusters for k-means
  double threshold{4.0};        // Threshold for hierarchical clustering
  int connectivity{6};          // 6 or 26 connectivity
  int minClusterSize{5};        // Minimum cluster size for HDBSCAN
  int minSamples{5};            // Minimum samples for HDBSCAN core points
};

// Run a clustering algorithm on the given occupancy grid
// Returns the clustering result
// The algorithms operate directly on the grid data to keep implementations pure
ClusterResult runClusterAlgorithm(
    ClusterAlgorithm algo,
    const ClusterParams& params,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// Individual algorithm implementations (kept simple for benchmarking fairness)

// Traditional DFS clustering - O(n) where n is number of occupied voxels
ClusterResult traditionalDFS(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// Skip-DFS clustering - O(n/skip^3) expected
ClusterResult skipDFS(
    int nx, int ny, int nz, int skip,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// DBSCAN clustering with grid-based spatial indexing
ClusterResult dbscan(
    int nx, int ny, int nz,
    double eps, int minPts,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// Hierarchical (single-linkage) clustering
ClusterResult hierarchical(
    int nx, int ny, int nz,
    double threshold,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// K-means clustering
ClusterResult kmeans(
    int nx, int ny, int nz,
    int k,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// Simplified spectral clustering (without eigendecomposition library)
ClusterResult spectral(
    int nx, int ny, int nz,
    int k,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// GCBD (Grid-based Connectivity using Union-Find)
ClusterResult gcbd(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// HDBSCAN (Hierarchical Density-Based Spatial Clustering)
ClusterResult hdbscan(
    int nx, int ny, int nz,
    int minClusterSize, int minSamples,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// CC3D (Connected Components 3D) - Fair basic implementation
// Uses basic Union-Find WITHOUT path compression or union-by-rank
// This ensures fair comparison with BLS at equivalent optimization levels
ClusterResult cc3d(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

// CC3D Optimized - Uses path compression and union-by-rank
// For reference comparison only (not for fair benchmarking against BLS)
ClusterResult cc3dOptimized(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

}  // namespace bls
