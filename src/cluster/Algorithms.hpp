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
  GCBD,          // union-find CCL; WITHDRAWN from the comparison (author 27.09.26): not Du & Wu's GCBD
  HDBSCAN,       // Hierarchical DBSCAN
  CC3D,          // Connected Components 3D (fair: basic Union-Find)
  // Dual-track pairs. Each <X> is the textbook-fair variant held to the same
  // optimization level as BLS; each <X>Optimized is the method as it is
  // actually published and deployed. Both tracks are reported.
  CC3DOptimized, // CC3D as SAUF (Wu, Otoo, Shoshani 2005; Wu, Otoo, Suzuki 2009)
  RLECCL,        // Run-Length Encoding CCL (textbook: per-voxel union-find)
  RLECCLOptimized,  // RLE-CCL as published (He, Chao, Suzuki 2008)
  VCCS,          // Voxel Cloud Connected Segmentation (textbook: uniform seeds)
  VCCSOptimized  // VCCS as published (Papon et al., CVPR 2013)
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

  // Seeding instrumentation. Written only by the two VCCS tracks and left at
  // zero by every other algorithm; nothing branches on them. They exist so the
  // seed budget can be audited against the code that actually runs, rather
  // than by reimplementing the seeding rules somewhere they can drift.
  //   seedCandidates     - fair: seed-grid points examined.
  //                        optimized: seed cells containing structure.
  //   seedsPlaced        - candidates that became a supervoxel seed.
  //   seedPruneThreshold - optimized: the paper's filter, pi*R_search^2 voxels
  //                        (a planar slice through the search volume, R_search =
  //                        S/2); a candidate survives iff at least this many
  //                        occupied voxels lie within R_search. 0 in the fair
  //                        track, which does not prune.
  int seedCandidates{0};
  int seedsPlaced{0};
  double seedPruneThreshold{0.0};
};

// Parameters for clustering algorithms
struct ClusterParams {
  int nx{0}, ny{0}, nz{0};     // Grid dimensions
  // Named for what it feeds, after a single `skip` fed two unrelated
  // algorithms and gave the (since removed) octree_ccl a leaf size of 3
  // against its documented intent of 8.
  int skipDfsJumpDistance{3};   // Jump distance for cluster::skipDFS
  double eps{3.0};              // Epsilon for DBSCAN
  int minPts{10};               // MinPts for DBSCAN
  int k{20};                    // Number of clusters for k-means
  double threshold{4.0};        // Threshold for hierarchical clustering
  int connectivity{6};          // 6 or 26 connectivity
  int minClusterSize{5};        // Minimum cluster size for HDBSCAN
  int minSamples{5};            // Minimum samples for HDBSCAN core points
  // Deck PBC xyz: adjacency wraps across the grid faces. Only methods for which
  // supportsPeriodic() is true accept it; runClusterAlgorithm throws for the others.
  bool periodic{false};
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
// CC3DOptimized, GCBD, RLECCL -- and, through Analyzer, by BLS.
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

// Whether the method implements periodic adjacency (deck PBC xyz).
bool supportsPeriodic(ClusterAlgorithm algo);

// Individual algorithm implementations (kept simple for benchmarking fairness)

// Traditional DFS clustering - O(n) where n is number of occupied voxels
ClusterResult traditionalDFS(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr,
    bool periodic = false);

// Skip-DFS clustering - O(n/skip^3) expected
ClusterResult skipDFS(
    int nx, int ny, int nz, int skip,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr,
    bool periodic = false);

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
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// K-means clustering
ClusterResult kmeans(
    int nx, int ny, int nz,
    int k,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited);


// GCBD (Grid-based Connectivity using Union-Find) -- WITHDRAWN from the comparison (author,
// 27.09.26): a plain 6-connected union-find CCL, identical in output to CC3D, and not the
// density-based grid clustering of Du & Wu (Entropy 24:1606, 2022) it was attributed to.
// Kept for the record and the exact-partition tests; not run in the campaigns.
ClusterResult gcbd(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr,
    bool periodic = false);

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

// CC3D Optimized - SAUF (Scan plus Array-based Union-Find) as published: Wu, Otoo &
// Shoshani, SPIE 5747 (2005); Wu, Otoo & Suzuki, Pattern Anal. Appl. 12:117-135 (2009).
// Scanning phase (S0, eq. 8) with the papers' array union-find, flattenL analysis phase,
// labeling pass over the whole label image. See the definition.
ClusterResult cc3dOptimized(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr,
    bool periodic = false);

// RLE-based CCL, textbook track. Encodes runs but unions individual voxels, and
// sweeps the full grid volume three times. Held at BLS's optimization level.
ClusterResult rleCCL(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr);

// RLE-CCL Optimized - He, Chao & Suzuki, "A run-based two-scan labeling algorithm",
// IEEE Trans. Image Process. 17(5):749-756, 2008, as published, in 3-D: first scan
// labelling runs with provisional label lists (rtable, next, tail), second scan over the
// whole label image. See the definition.
ClusterResult rleCCLOptimized(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels = nullptr,
    bool periodic = false);

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

// VCCS Optimized - Voxel Cloud Connectivity Segmentation as published: Papon, Abramov,
// Schoeler & Woergoetter, CVPR 2013, doi:10.1109/CVPR.2013.264, §3.1-3.4 (audit decision D6
// revised): 26-adjacency, nearest-to-centre seeding, planar-slice seed filter, BFS flow
// with per-level interleaving of supervoxels, centre updates, at most 5 iterations. Spatial
// distance term only (no colour or normals on an occupancy grid). See the definition.
ClusterResult vccsOptimized(
    int nx, int ny, int nz,
    double seedResolution,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    bool periodic = false);

}  // namespace bls
