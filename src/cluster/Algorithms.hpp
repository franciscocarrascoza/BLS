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
  CC3DOptimized, // Connected Components 3D (optimized: path compression + union-by-rank)
  RLECCL,        // Run-Length Encoding CCL (sparse-aware, Category 2a)
  OctreeCCL,     // Octree-based CCL (hierarchical, skips empty octants, Category 2a)
  VCCS           // Voxel Cloud Connected Segmentation (uniform-grid seeding, Category 2e)
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

// Optional per-voxel label output.
//
// Comparing cluster counts alone cannot distinguish a correct labelling from
// one that merges one pair of components and splits another; the counts match
// while the partition is wrong. Differential testing therefore needs the
// partition itself, not a summary of it.
//
// Contract when a caller passes a non-null `labels`:
//   - it is resized to nx*ny*nz and overwritten in full;
//   - unoccupied voxels get -1;
//   - occupied voxels get dense ids 0..nclusters-1.
// The particular id an algorithm assigns to a given component is arbitrary and
// differs between algorithms. Canonical relabelling for comparison is the
// harness's job, deliberately not done here: forcing a canonical order inside
// each algorithm would cost a sort that the timed benchmark path must not pay,
// and would hide exactly the ordering differences a harness may want to see.
//
// Passing nullptr (the default) allocates nothing and leaves the benchmark
// path as it was.
//
// Supported by the six exact algorithms -- TraditionalDFS, CC3D,
// CC3DOptimized, GCBD, RLECCL, OctreeCCL -- and, through Analyzer, by BLS.
// Requesting labels from any other algorithm throws rather than returning a
// buffer that silently does not mean what it appears to.

// Run a clustering algorithm on the given occupancy grid
// Returns the clustering result
// The algorithms operate directly on the grid data to keep implementations pure
ClusterResult runClusterAlgorithm(
    ClusterAlgorithm algo,
    const ClusterParams& params,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// True if runClusterAlgorithm accepts a non-null `labels` for this algorithm.
bool supportsLabels(ClusterAlgorithm algo);

// Individual algorithm implementations (kept simple for benchmarking fairness)

// Traditional DFS clustering - O(n) where n is number of occupied voxels
ClusterResult traditionalDFS(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

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
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

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
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// CC3D Optimized - Uses path compression and union-by-rank
// For reference comparison only (not for fair benchmarking against BLS)
ClusterResult cc3dOptimized(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// RLE-based CCL - Run-Length Encoding Connected Component Labeling
// Encodes consecutive occupied voxels as runs; merges adjacent runs via union-find.
// Sparse-aware: complexity scales with number of runs, not grid volume.
ClusterResult rleCCL(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// Octree-based CCL - Hierarchical grid subdivision
// Recursively skips empty octants; runs local connectivity at leaves.
// Shares BLS's "skip empty space" philosophy via a different strategy.
// leafSize: stop subdividing when all dimensions <= leafSize (default: 8)
ClusterResult octreeCCL(
    int nx, int ny, int nz,
    int leafSize,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// VCCS - Voxel Cloud Connected Segmentation
// Places seeds on a uniform 3D grid (spacing = seedResolution voxels) then
// expands regions outward via BFS ordered by distance to seed.
// Structurally closest to BLS (seed + expand), but uses uniform rather than
// crystallographic seeding — does NOT guarantee topological correctness.
ClusterResult vccs(
    int nx, int ny, int nz,
    double seedResolution,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);

}  // namespace bls
