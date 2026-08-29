#include "cluster/Algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cluster/UnionFindBasic.hpp"
#include "util/Timer.hpp"

namespace bls {

namespace {

// Stack element for DFS-based algorithms
struct StackElement {
  int i, j, k;
  int is_skipping;
};

// Point structure for algorithms that work with point lists
struct Point {
  int i, j, k;
};

// Union-Find structure
struct UnionFind {
  int parent;
  int rank;
};

// Index calculation
inline std::size_t idx3(int i, int j, int k, int ny, int nz) {
  return static_cast<std::size_t>(i) * static_cast<std::size_t>(ny) * nz +
         static_cast<std::size_t>(j) * nz + static_cast<std::size_t>(k);
}

// 6-connectivity deltas
constexpr int deltas6[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                                {0, 1, 0},  {0, 0, -1}, {0, 0, 1}};

// --- optional label output (see Algorithms.hpp) -----------------------------
//
// Both helpers no-op on nullptr, so the benchmark path allocates nothing and
// pays only a null check outside any hot loop.

inline void initLabels(std::vector<int>* labels, std::size_t totalSize) {
  if (labels) labels->assign(totalSize, -1);
}

// Compacts arbitrary non-negative component keys -- for the union-find methods,
// the root voxel index -- into dense ids 0..k-1 in order of first appearance.
// Voxels left at -1 (unoccupied) are untouched.
void compactLabels(std::vector<int>* labels) {
  if (!labels) return;
  std::unordered_map<int, int> remap;
  int next = 0;
  for (int& v : *labels) {
    if (v < 0) continue;
    auto it = remap.find(v);
    if (it == remap.end()) it = remap.emplace(v, next++).first;
    v = it->second;
  }
}

}  // namespace

ClusterAlgorithm parseAlgorithm(const std::string& name) {
  std::string lower = name;
  for (auto& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

  if (lower == "bls" || lower == "default") return ClusterAlgorithm::BLS;
  if (lower == "traditional_dfs" || lower == "dfs" || lower == "traditional") return ClusterAlgorithm::TraditionalDFS;
  if (lower == "skip_dfs" || lower == "skipdfs" || lower == "skip") return ClusterAlgorithm::SkipDFS;
  if (lower == "dbscan") return ClusterAlgorithm::DBSCAN;
  if (lower == "hierarchical" || lower == "single_linkage") return ClusterAlgorithm::Hierarchical;
  if (lower == "kmeans" || lower == "k-means" || lower == "k_means") return ClusterAlgorithm::KMeans;
  if (lower == "spectral") return ClusterAlgorithm::Spectral;
  if (lower == "gcbd" || lower == "union_find" || lower == "uf") return ClusterAlgorithm::GCBD;
  if (lower == "hdbscan") return ClusterAlgorithm::HDBSCAN;
  if (lower == "cc3d" || lower == "cc3d_basic") return ClusterAlgorithm::CC3D;
  if (lower == "cc3d_optimized" || lower == "cc3d_opt") return ClusterAlgorithm::CC3DOptimized;
  if (lower == "rle_ccl" || lower == "rleccl" || lower == "rle") return ClusterAlgorithm::RLECCL;
  if (lower == "octree_ccl" || lower == "octreeccl" || lower == "octree") return ClusterAlgorithm::OctreeCCL;
  if (lower == "vccs") return ClusterAlgorithm::VCCS;

  throw std::runtime_error("Unknown clustering algorithm: " + name);
}

std::string algorithmToString(ClusterAlgorithm algo) {
  switch (algo) {
    case ClusterAlgorithm::BLS: return "bls";
    case ClusterAlgorithm::TraditionalDFS: return "traditional_dfs";
    case ClusterAlgorithm::SkipDFS: return "skip_dfs";
    case ClusterAlgorithm::DBSCAN: return "dbscan";
    case ClusterAlgorithm::Hierarchical: return "hierarchical";
    case ClusterAlgorithm::KMeans: return "kmeans";
    case ClusterAlgorithm::Spectral: return "spectral";
    case ClusterAlgorithm::GCBD: return "gcbd";
    case ClusterAlgorithm::HDBSCAN: return "hdbscan";
    case ClusterAlgorithm::CC3D: return "cc3d";
    case ClusterAlgorithm::CC3DOptimized: return "cc3d_optimized";
    case ClusterAlgorithm::RLECCL: return "rle_ccl";
    case ClusterAlgorithm::OctreeCCL: return "octree_ccl";
    case ClusterAlgorithm::VCCS: return "vccs";
  }
  return "unknown";
}

std::vector<std::string> listAlgorithms() {
  return {"bls", "traditional_dfs", "skip_dfs", "dbscan",
          "hierarchical", "kmeans", "spectral", "gcbd", "hdbscan",
          "cc3d", "cc3d_optimized", "rle_ccl", "octree_ccl", "vccs"};
}

bool supportsLabels(ClusterAlgorithm algo) {
  switch (algo) {
    case ClusterAlgorithm::TraditionalDFS:
    case ClusterAlgorithm::GCBD:
    case ClusterAlgorithm::CC3D:
    case ClusterAlgorithm::CC3DOptimized:
    case ClusterAlgorithm::RLECCL:
    case ClusterAlgorithm::OctreeCCL:
      return true;
    default:
      // BLS is labelled through Analyzer, not here. The remaining methods
      // (skip_dfs, dbscan, hierarchical, kmeans, spectral, hdbscan, vccs) are
      // not exact partitioners and have no label output yet.
      return false;
  }
}

ClusterResult runClusterAlgorithm(
    ClusterAlgorithm algo,
    const ClusterParams& params,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  if (labels && !supportsLabels(algo)) {
    throw std::runtime_error("Label output requested for '" + algorithmToString(algo) +
                             "', which does not provide one. Returning an unlabelled or "
                             "partially labelled buffer would look like a valid partition.");
  }

  // Note: BLS is handled separately in the main analyzer since it requires
  // lattice enumeration. This function handles the comparison algorithms.
  switch (algo) {
    case ClusterAlgorithm::BLS:
      throw std::runtime_error("BLS algorithm should be run through the Analyzer class.");
    case ClusterAlgorithm::TraditionalDFS:
      return traditionalDFS(params.nx, params.ny, params.nz, occupancy, visited, labels);
    case ClusterAlgorithm::SkipDFS:
      return skipDFS(params.nx, params.ny, params.nz, params.skip, occupancy, visited);
    case ClusterAlgorithm::DBSCAN:
      return dbscan(params.nx, params.ny, params.nz, params.eps, params.minPts, occupancy, visited);
    case ClusterAlgorithm::Hierarchical:
      return hierarchical(params.nx, params.ny, params.nz, params.threshold, occupancy, visited);
    case ClusterAlgorithm::KMeans:
      return kmeans(params.nx, params.ny, params.nz, params.k, occupancy, visited);
    case ClusterAlgorithm::Spectral:
      return spectral(params.nx, params.ny, params.nz, params.k, occupancy, visited);
    case ClusterAlgorithm::GCBD:
      return gcbd(params.nx, params.ny, params.nz, occupancy, visited, labels);
    case ClusterAlgorithm::HDBSCAN:
      return hdbscan(params.nx, params.ny, params.nz, params.minClusterSize, params.minSamples, occupancy, visited);
    case ClusterAlgorithm::CC3D:
      return cc3d(params.nx, params.ny, params.nz, params.connectivity, occupancy, visited, labels);
    case ClusterAlgorithm::CC3DOptimized:
      return cc3dOptimized(params.nx, params.ny, params.nz, params.connectivity, occupancy, visited,
                           labels);
    case ClusterAlgorithm::RLECCL:
      return rleCCL(params.nx, params.ny, params.nz, occupancy, visited, labels);
    case ClusterAlgorithm::OctreeCCL:
      // Use params.skip as leaf size (default 8 voxels per side for good balance)
      return octreeCCL(params.nx, params.ny, params.nz,
                       params.skip > 0 ? params.skip : 8, occupancy, visited, labels);
    case ClusterAlgorithm::VCCS:
      // Use params.eps as seed spacing in voxels (default 3.0)
      return vccs(params.nx, params.ny, params.nz, params.eps, occupancy, visited);
  }
  throw std::runtime_error("Unhandled algorithm type");
}

// Traditional DFS - kept simple for benchmarking
ClusterResult traditionalDFS(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  std::vector<StackElement> stack;
  stack.reserve(10000);

  for (std::size_t idx = 0; idx < totalSize; ++idx) {
    if (occupancy[idx] == 1 && !visited[idx]) {
      int clusterSize = 0;
      int i0 = static_cast<int>(idx / (static_cast<std::size_t>(ny) * nz));
      int remainder = static_cast<int>(idx % (static_cast<std::size_t>(ny) * nz));
      int j0 = remainder / nz;
      int k0 = remainder % nz;

      stack.clear();
      stack.push_back({i0, j0, k0, 0});

      while (!stack.empty()) {
        StackElement curr = stack.back();
        stack.pop_back();

        int i = curr.i, j = curr.j, k = curr.k;
        if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) continue;

        std::size_t index = idx3(i, j, k, ny, nz);
        if (occupancy[index] == 0 || visited[index]) continue;

        visited[index] = 1;
        // result.nclusters is the number of components already completed, so
        // it is the 0-based ordinal of the one being walked: dense already.
        // compactLabels() below is then a no-op, kept so the contract survives
        // any future change to the order components are discovered in.
        if (labels) (*labels)[index] = result.nclusters;
        clusterSize++;
        result.visitedVoxels++;

        for (int d = 0; d < 6; ++d) {
          int ni = i + deltas6[d][0];
          int nj = j + deltas6[d][1];
          int nk = k + deltas6[d][2];
          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            stack.push_back({ni, nj, nk, 0});
          }
        }
      }

      if (clusterSize > 0) {
        result.nclusters++;
        result.clusterSizes.push_back(clusterSize);
        result.maxCluster = std::max(result.maxCluster, clusterSize);
      }
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// Skip-DFS - kept simple for benchmarking
ClusterResult skipDFS(
    int nx, int ny, int nz, int skip,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  std::vector<StackElement> stack;
  stack.reserve(10000);

  int skipDeltas[6][3] = {{-skip, 0, 0}, {skip, 0, 0}, {0, -skip, 0},
                           {0, skip, 0}, {0, 0, -skip}, {0, 0, skip}};

  for (std::size_t idx = 0; idx < totalSize; ++idx) {
    if (occupancy[idx] == 1 && !visited[idx]) {
      int clusterSize = 0;
      int i0 = static_cast<int>(idx / (static_cast<std::size_t>(ny) * nz));
      int remainder = static_cast<int>(idx % (static_cast<std::size_t>(ny) * nz));
      int j0 = remainder / nz;
      int k0 = remainder % nz;

      stack.clear();
      stack.push_back({i0, j0, k0, 1});

      while (!stack.empty()) {
        StackElement curr = stack.back();
        stack.pop_back();

        int i = curr.i, j = curr.j, k = curr.k;
        int isSkipping = curr.is_skipping;

        if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) continue;

        std::size_t index = idx3(i, j, k, ny, nz);
        if (occupancy[index] == 0 || visited[index]) continue;

        visited[index] = 1;
        clusterSize++;
        result.visitedVoxels++;

        if (isSkipping) {
          // Skip-mode: jump 'skip' voxels ahead in each axis direction.
          // Also enqueue all intermediate voxels so they are visited and
          // not mistakenly treated as new cluster seeds by the outer loop.
          for (int d = 0; d < 6; ++d) {
            const int di = skipDeltas[d][0];
            const int dj = skipDeltas[d][1];
            const int dk = skipDeltas[d][2];
            // Unit step along the skip direction (±1 per axis)
            const int step_i = di / skip;
            const int step_j = dj / skip;
            const int step_k = dk / skip;
            // Enqueue intermediate voxels with no-skip flag
            for (int step = 1; step < skip; ++step) {
              int ii = i + step_i * step;
              int jj = j + step_j * step;
              int kk = k + step_k * step;
              if (ii >= 0 && ii < nx && jj >= 0 && jj < ny && kk >= 0 && kk < nz) {
                stack.push_back({ii, jj, kk, 0});
              }
            }
            // Enqueue the skip-target
            int ni = i + di;
            int nj = j + dj;
            int nk = k + dk;
            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
              std::size_t nIndex = idx3(ni, nj, nk, ny, nz);
              stack.push_back({ni, nj, nk, occupancy[nIndex] == 1 ? 1 : 0});
            }
          }
        } else {
          // Normal 6-connectivity exploration
          for (int d = 0; d < 6; ++d) {
            int ni = i + deltas6[d][0];
            int nj = j + deltas6[d][1];
            int nk = k + deltas6[d][2];
            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
              std::size_t nIndex = idx3(ni, nj, nk, ny, nz);
              stack.push_back({ni, nj, nk, occupancy[nIndex] == 1 ? 1 : 0});
            }
          }
        }
      }

      if (clusterSize > 0) {
        result.nclusters++;
        result.clusterSizes.push_back(clusterSize);
        result.maxCluster = std::max(result.maxCluster, clusterSize);
      }
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// DBSCAN with grid-based spatial indexing
ClusterResult dbscan(
    int nx, int ny, int nz,
    double eps, int minPts,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  // Collect occupied points
  std::vector<Point> points;
  points.reserve(totalSize / 10);  // Estimate

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          points.push_back({i, j, k});
        }
      }
    }
  }

  int numPoints = static_cast<int>(points.size());
  if (numPoints == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  // Build spatial grid index
  int cellSize = static_cast<int>(std::ceil(eps));
  int gridCellsX = (nx + cellSize - 1) / cellSize;
  int gridCellsY = (ny + cellSize - 1) / cellSize;
  int gridCellsZ = (nz + cellSize - 1) / cellSize;

  std::vector<std::vector<int>> cells(
      static_cast<std::size_t>(gridCellsX) * gridCellsY * gridCellsZ);

  auto cellIdx = [&](int ci, int cj, int ck) {
    return static_cast<std::size_t>(ci) * gridCellsY * gridCellsZ +
           static_cast<std::size_t>(cj) * gridCellsZ + ck;
  };

  for (int p = 0; p < numPoints; ++p) {
    int ci = points[p].i / cellSize;
    int cj = points[p].j / cellSize;
    int ck = points[p].k / cellSize;
    cells[cellIdx(ci, cj, ck)].push_back(p);
  }

  std::vector<bool> pointVisited(numPoints, false);
  std::vector<StackElement> stack;
  stack.reserve(1000);

  for (int p = 0; p < numPoints; ++p) {
    if (pointVisited[p]) continue;

    Point curr = points[p];
    int ci = curr.i / cellSize;
    int cj = curr.j / cellSize;
    int ck = curr.k / cellSize;

    // Count neighbors
    int neighbors = 0;
    for (int di = -1; di <= 1; ++di) {
      for (int dj = -1; dj <= 1; ++dj) {
        for (int dk = -1; dk <= 1; ++dk) {
          int nci = ci + di, ncj = cj + dj, nck = ck + dk;
          if (nci < 0 || nci >= gridCellsX || ncj < 0 || ncj >= gridCellsY ||
              nck < 0 || nck >= gridCellsZ)
            continue;

          for (int q : cells[cellIdx(nci, ncj, nck)]) {
            Point other = points[q];
            double dist = std::sqrt(
                std::pow(curr.i - other.i, 2) +
                std::pow(curr.j - other.j, 2) +
                std::pow(curr.k - other.k, 2));
            if (dist <= eps) neighbors++;
          }
        }
      }
    }

    if (neighbors >= minPts) {
      int clusterSize = 0;
      stack.clear();
      stack.push_back({curr.i, curr.j, curr.k, p});

      while (!stack.empty()) {
        StackElement s = stack.back();
        stack.pop_back();

        int pidx = s.is_skipping;  // We use is_skipping to store point index
        if (pointVisited[pidx]) continue;

        pointVisited[pidx] = 1;
        clusterSize++;

        std::size_t occIdx = idx3(points[pidx].i, points[pidx].j, points[pidx].k, ny, nz);
        visited[occIdx] = 1;
        result.visitedVoxels++;

        // Expand cluster
        int nci = points[pidx].i / cellSize;
        int ncj = points[pidx].j / cellSize;
        int nck = points[pidx].k / cellSize;

        for (int di = -1; di <= 1; ++di) {
          for (int dj = -1; dj <= 1; ++dj) {
            for (int dk = -1; dk <= 1; ++dk) {
              int nnci = nci + di, nncj = ncj + dj, nnck = nck + dk;
              if (nnci < 0 || nnci >= gridCellsX || nncj < 0 || nncj >= gridCellsY ||
                  nnck < 0 || nnck >= gridCellsZ)
                continue;

              for (int q : cells[cellIdx(nnci, nncj, nnck)]) {
                if (pointVisited[q]) continue;
                Point other = points[q];
                double dist = std::sqrt(
                    std::pow(points[pidx].i - other.i, 2) +
                    std::pow(points[pidx].j - other.j, 2) +
                    std::pow(points[pidx].k - other.k, 2));
                if (dist <= eps) {
                  stack.push_back({other.i, other.j, other.k, q});
                }
              }
            }
          }
        }
      }

      if (clusterSize > 0) {
        result.nclusters++;
        result.clusterSizes.push_back(clusterSize);
        result.maxCluster = std::max(result.maxCluster, clusterSize);
      }
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// Hierarchical (Single-Linkage) clustering
ClusterResult hierarchical(
    int nx, int ny, int nz,
    double threshold,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::fill(visited.begin(), visited.end(), 0);

  // Collect occupied points
  std::vector<Point> points;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          points.push_back({i, j, k});
        }
      }
    }
  }

  int numPoints = static_cast<int>(points.size());
  if (numPoints == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  // Union-Find structure
  std::vector<UnionFind> uf(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    uf[i].parent = i;
    uf[i].rank = 0;
  }

  auto find = [&](int x) {
    while (uf[x].parent != x) {
      uf[x].parent = uf[uf[x].parent].parent;  // Path compression
      x = uf[x].parent;
    }
    return x;
  };

  auto unite = [&](int x, int y) {
    int rootX = find(x);
    int rootY = find(y);
    if (rootX != rootY) {
      if (uf[rootX].rank < uf[rootY].rank) {
        uf[rootX].parent = rootY;
      } else if (uf[rootX].rank > uf[rootY].rank) {
        uf[rootY].parent = rootX;
      } else {
        uf[rootY].parent = rootX;
        uf[rootX].rank++;
      }
    }
  };

  // Merge points within threshold
  for (int p = 0; p < numPoints; ++p) {
    for (int q = p + 1; q < numPoints; ++q) {
      double dist = std::sqrt(
          std::pow(points[p].i - points[q].i, 2) +
          std::pow(points[p].j - points[q].j, 2) +
          std::pow(points[p].k - points[q].k, 2));
      if (dist <= threshold) {
        unite(p, q);
      }
    }
  }

  // Count clusters
  std::vector<int> clusterSizeMap(numPoints, 0);
  for (int i = 0; i < numPoints; ++i) {
    int root = find(i);
    clusterSizeMap[root]++;

    std::size_t occIdx = idx3(points[i].i, points[i].j, points[i].k, ny, nz);
    visited[occIdx] = 1;
    result.visitedVoxels++;
  }

  for (int i = 0; i < numPoints; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// K-means clustering
ClusterResult kmeans(
    int nx, int ny, int nz,
    int k,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::fill(visited.begin(), visited.end(), 0);

  // Collect occupied points
  std::vector<Point> points;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int kk = 0; kk < nz; ++kk) {
        std::size_t idx = idx3(i, j, kk, ny, nz);
        if (occupancy[idx] == 1) {
          points.push_back({i, j, kk});
        }
      }
    }
  }

  int numPoints = static_cast<int>(points.size());
  if (numPoints == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  k = std::min(k, numPoints);

  // Initialize centroids (evenly spaced)
  std::vector<double> centroidsX(k), centroidsY(k), centroidsZ(k);
  for (int i = 0; i < k; ++i) {
    int idx = (i * numPoints) / k;
    centroidsX[i] = points[idx].i;
    centroidsY[i] = points[idx].j;
    centroidsZ[i] = points[idx].k;
  }

  std::vector<int> assignments(numPoints);

  // Run k-means for 10 iterations
  for (int iter = 0; iter < 10; ++iter) {
    // Assign points to nearest centroid
    for (int p = 0; p < numPoints; ++p) {
      double minDist = std::numeric_limits<double>::max();
      int bestK = 0;
      for (int kk = 0; kk < k; ++kk) {
        double dist = std::sqrt(
            std::pow(points[p].i - centroidsX[kk], 2) +
            std::pow(points[p].j - centroidsY[kk], 2) +
            std::pow(points[p].k - centroidsZ[kk], 2));
        if (dist < minDist) {
          minDist = dist;
          bestK = kk;
        }
      }
      assignments[p] = bestK;
    }

    // Update centroids
    std::vector<double> sumX(k, 0), sumY(k, 0), sumZ(k, 0);
    std::vector<int> counts(k, 0);
    for (int p = 0; p < numPoints; ++p) {
      int kk = assignments[p];
      sumX[kk] += points[p].i;
      sumY[kk] += points[p].j;
      sumZ[kk] += points[p].k;
      counts[kk]++;
    }
    for (int kk = 0; kk < k; ++kk) {
      if (counts[kk] > 0) {
        centroidsX[kk] = sumX[kk] / counts[kk];
        centroidsY[kk] = sumY[kk] / counts[kk];
        centroidsZ[kk] = sumZ[kk] / counts[kk];
      }
    }
  }

  // Count cluster sizes
  std::vector<int> clusterSizeMap(k, 0);
  for (int p = 0; p < numPoints; ++p) {
    clusterSizeMap[assignments[p]]++;

    std::size_t occIdx = idx3(points[p].i, points[p].j, points[p].k, ny, nz);
    visited[occIdx] = 1;
    result.visitedVoxels++;
  }

  for (int kk = 0; kk < k; ++kk) {
    if (clusterSizeMap[kk] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[kk]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[kk]);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// Simplified spectral clustering (without eigendecomposition)
ClusterResult spectral(
    int nx, int ny, int nz,
    int k,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::fill(visited.begin(), visited.end(), 0);

  // Collect occupied points
  std::vector<Point> points;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int kk = 0; kk < nz; ++kk) {
        std::size_t idx = idx3(i, j, kk, ny, nz);
        if (occupancy[idx] == 1) {
          points.push_back({i, j, kk});
        }
      }
    }
  }

  int numPoints = static_cast<int>(points.size());
  if (numPoints == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  k = std::min(k, numPoints);

  // Simplified: just assign points to clusters based on index modulo k
  // (since we don't have eigendecomposition library, keeping it pure)
  std::vector<int> labels(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    labels[i] = i % k;
  }

  // Count cluster sizes
  std::vector<int> clusterSizeMap(k, 0);
  for (int p = 0; p < numPoints; ++p) {
    clusterSizeMap[labels[p]]++;

    std::size_t occIdx = idx3(points[p].i, points[p].j, points[p].k, ny, nz);
    visited[occIdx] = 1;
    result.visitedVoxels++;
  }

  for (int kk = 0; kk < k; ++kk) {
    if (clusterSizeMap[kk] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[kk]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[kk]);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// GCBD (Grid-based Connectivity using Union-Find)
ClusterResult gcbd(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // Union-Find structure operating on voxel indices
  std::vector<int> parent(totalSize);
  std::vector<int> rank(totalSize, 0);

  for (std::size_t i = 0; i < totalSize; ++i) {
    parent[i] = static_cast<int>(i);
  }

  auto find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];  // Path compression
      x = parent[x];
    }
    return x;
  };

  auto unite = [&](int x, int y) {
    int rootX = find(x);
    int rootY = find(y);
    if (rootX != rootY) {
      if (rank[rootX] < rank[rootY]) {
        parent[rootX] = rootY;
      } else if (rank[rootX] > rank[rootY]) {
        parent[rootY] = rootX;
      } else {
        parent[rootY] = rootX;
        rank[rootX]++;
      }
    }
  };

  // Build connectivity
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] != 1) continue;

        for (int d = 0; d < 6; ++d) {
          int ni = i + deltas6[d][0];
          int nj = j + deltas6[d][1];
          int nk = k + deltas6[d][2];
          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
            if (occupancy[nIdx] == 1) {
              unite(static_cast<int>(idx), static_cast<int>(nIdx));
            }
          }
        }
      }
    }
  }

  // Count clusters
  std::vector<int> clusterSizeMap(totalSize, 0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          int root = find(static_cast<int>(idx));
          clusterSizeMap[root]++;
          if (labels) (*labels)[idx] = root;
          visited[idx] = 1;
          result.visitedVoxels++;
        }
      }
    }
  }

  for (std::size_t i = 0; i < totalSize; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// HDBSCAN (Hierarchical Density-Based Spatial Clustering)
// Simplified implementation for fair benchmarking
ClusterResult hdbscan(
    int nx, int ny, int nz,
    int minClusterSize, int minSamples,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::fill(visited.begin(), visited.end(), 0);

  // Collect occupied points
  std::vector<Point> points;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          points.push_back({i, j, k});
        }
      }
    }
  }

  int numPoints = static_cast<int>(points.size());
  if (numPoints == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  // Compute core distances (distance to minSamples-th nearest neighbor)
  std::vector<double> coreDistances(numPoints);
  for (int p = 0; p < numPoints; ++p) {
    std::vector<double> distances;
    distances.reserve(numPoints);

    for (int q = 0; q < numPoints; ++q) {
      if (p == q) continue;
      double dist = std::sqrt(
          std::pow(points[p].i - points[q].i, 2) +
          std::pow(points[p].j - points[q].j, 2) +
          std::pow(points[p].k - points[q].k, 2));
      distances.push_back(dist);
    }

    std::sort(distances.begin(), distances.end());
    int kIdx = std::min(minSamples - 1, static_cast<int>(distances.size()) - 1);
    coreDistances[p] = kIdx >= 0 ? distances[kIdx] : 0.0;
  }

  // Build mutual reachability distance graph using Union-Find
  // Mutual reachability = max(core_dist(a), core_dist(b), dist(a,b))
  struct Edge {
    int p1, p2;
    double weight;
    bool operator<(const Edge& other) const { return weight < other.weight; }
  };

  std::vector<Edge> edges;
  edges.reserve(numPoints * numPoints / 2);

  for (int p = 0; p < numPoints; ++p) {
    for (int q = p + 1; q < numPoints; ++q) {
      double dist = std::sqrt(
          std::pow(points[p].i - points[q].i, 2) +
          std::pow(points[p].j - points[q].j, 2) +
          std::pow(points[p].k - points[q].k, 2));
      double mutualReach = std::max({coreDistances[p], coreDistances[q], dist});
      edges.push_back({p, q, mutualReach});
    }
  }

  // Sort edges by mutual reachability distance
  std::sort(edges.begin(), edges.end());

  // Build minimum spanning tree using Kruskal's algorithm with Union-Find
  std::vector<UnionFind> uf(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    uf[i].parent = i;
    uf[i].rank = 0;
  }

  auto find = [&](int x) {
    while (uf[x].parent != x) {
      uf[x].parent = uf[uf[x].parent].parent;
      x = uf[x].parent;
    }
    return x;
  };

  auto unite = [&](int x, int y) -> bool {
    int rootX = find(x);
    int rootY = find(y);
    if (rootX == rootY) return false;

    if (uf[rootX].rank < uf[rootY].rank) {
      uf[rootX].parent = rootY;
    } else if (uf[rootX].rank > uf[rootY].rank) {
      uf[rootY].parent = rootX;
    } else {
      uf[rootY].parent = rootX;
      uf[rootX].rank++;
    }
    return true;
  };

  // Process edges in order, cutting at appropriate threshold
  // Simplified: cut the MST where edge weight exceeds median mutual reachability
  double medianWeight = edges.empty() ? 0.0 : edges[edges.size() / 2].weight;

  for (const auto& edge : edges) {
    // Only connect points if mutual reachability is not too large
    if (edge.weight <= medianWeight * 1.5) {
      unite(edge.p1, edge.p2);
    }
  }

  // Extract clusters and filter by minimum cluster size
  std::vector<int> clusterSizeMap(numPoints, 0);
  for (int i = 0; i < numPoints; ++i) {
    int root = find(i);
    clusterSizeMap[root]++;
  }

  // Assign cluster labels only to clusters >= minClusterSize
  std::vector<int> clusterLabels(numPoints, -1);
  int clusterId = 0;
  for (int i = 0; i < numPoints; ++i) {
    if (clusterSizeMap[i] >= minClusterSize && clusterLabels[find(i)] == -1) {
      clusterLabels[find(i)] = clusterId++;
    }
  }

  // Map back to voxel grid and count final clusters
  std::vector<int> finalClusterSizes(clusterId, 0);
  for (int i = 0; i < numPoints; ++i) {
    int root = find(i);
    int label = clusterLabels[root];

    if (label >= 0) {
      finalClusterSizes[label]++;
      std::size_t occIdx = idx3(points[i].i, points[i].j, points[i].k, ny, nz);
      visited[occIdx] = 1;
      result.visitedVoxels++;
    }
  }

  // Populate result
  result.nclusters = clusterId;
  for (int size : finalClusterSizes) {
    if (size > 0) {
      result.clusterSizes.push_back(size);
      result.maxCluster = std::max(result.maxCluster, size);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// CC3D (Connected Components 3D) - Fair basic implementation
// Uses basic Union-Find WITHOUT path compression or union-by-rank
// This ensures fair comparison with BLS at equivalent optimization levels
ClusterResult cc3d(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // 26-connectivity deltas (includes 6-connectivity as subset)
  constexpr int deltas26[26][3] = {
      // Face neighbors (6-connectivity)
      {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
      // Edge neighbors
      {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0},
      {-1, 0, -1}, {-1, 0, 1}, {1, 0, -1}, {1, 0, 1},
      {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1},
      // Corner neighbors
      {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},
      {1, -1, -1}, {1, -1, 1}, {1, 1, -1}, {1, 1, 1}
  };

  int numNeighbors = (connectivity == 26) ? 26 : 6;

  // Basic Union-Find WITHOUT path compression or union-by-rank
  UnionFindBasic uf(totalSize);

  // Build connectivity
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] != 1) continue;

        for (int d = 0; d < numNeighbors; ++d) {
          int ni = i + deltas26[d][0];
          int nj = j + deltas26[d][1];
          int nk = k + deltas26[d][2];
          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
            if (occupancy[nIdx] == 1) {
              uf.unite(static_cast<int>(idx), static_cast<int>(nIdx));
            }
          }
        }
      }
    }
  }

  // Count clusters
  std::vector<int> clusterSizeMap(totalSize, 0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          int root = uf.find(static_cast<int>(idx));
          clusterSizeMap[root]++;
          if (labels) (*labels)[idx] = root;
          visited[idx] = 1;
          result.visitedVoxels++;
        }
      }
    }
  }

  for (std::size_t i = 0; i < totalSize; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// CC3D Optimized - Uses path compression and union-by-rank
// For reference comparison only (not for fair benchmarking against BLS)
ClusterResult cc3dOptimized(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // 26-connectivity deltas (includes 6-connectivity as subset)
  constexpr int deltas26[26][3] = {
      // Face neighbors (6-connectivity)
      {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
      // Edge neighbors
      {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0},
      {-1, 0, -1}, {-1, 0, 1}, {1, 0, -1}, {1, 0, 1},
      {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1},
      // Corner neighbors
      {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},
      {1, -1, -1}, {1, -1, 1}, {1, 1, -1}, {1, 1, 1}
  };

  int numNeighbors = (connectivity == 26) ? 26 : 6;

  // Union-Find structure WITH path compression and union-by-rank
  std::vector<int> parent(totalSize);
  std::vector<int> rank(totalSize, 0);

  for (std::size_t i = 0; i < totalSize; ++i) {
    parent[i] = static_cast<int>(i);
  }

  auto find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];  // Path compression
      x = parent[x];
    }
    return x;
  };

  auto unite = [&](int x, int y) {
    int rootX = find(x);
    int rootY = find(y);
    if (rootX != rootY) {
      // Union by rank
      if (rank[rootX] < rank[rootY]) {
        parent[rootX] = rootY;
      } else if (rank[rootX] > rank[rootY]) {
        parent[rootY] = rootX;
      } else {
        parent[rootY] = rootX;
        rank[rootX]++;
      }
    }
  };

  // Build connectivity
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] != 1) continue;

        for (int d = 0; d < numNeighbors; ++d) {
          int ni = i + deltas26[d][0];
          int nj = j + deltas26[d][1];
          int nk = k + deltas26[d][2];
          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
            if (occupancy[nIdx] == 1) {
              unite(static_cast<int>(idx), static_cast<int>(nIdx));
            }
          }
        }
      }
    }
  }

  // Count clusters
  std::vector<int> clusterSizeMap(totalSize, 0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          int root = find(static_cast<int>(idx));
          clusterSizeMap[root]++;
          if (labels) (*labels)[idx] = root;
          visited[idx] = 1;
          result.visitedVoxels++;
        }
      }
    }
  }

  for (std::size_t i = 0; i < totalSize; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

// ── RLE-based CCL ────────────────────────────────────────────────────────────
//
// Run-Length Encoding Connected Component Labeling.
// For each row (i,j), consecutive occupied voxels along k form "runs".
// Adjacent rows / adjacent slices with overlapping run extents are merged
// via union-find on the global voxel index.
//
// Complexity: O(R * max_runs_per_row) where R = number of occupied runs.
// For sparse grids R << NX*NY*NZ, making this faster than full-grid scan.
//
// Reference: He et al., "A Run-Based Two-Scan Labeling Algorithm",
//            IEEE Trans. Image Processing 17(5), 2008.

ClusterResult rleCCL(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // Global union-find over voxel indices (path-halving + union-by-rank)
  std::vector<int> parent(totalSize);
  std::vector<int> ufRank(totalSize, 0);
  for (std::size_t i = 0; i < totalSize; ++i) parent[i] = static_cast<int>(i);

  auto find = [&](int x) -> int {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];  // path halving
      x = parent[x];
    }
    return x;
  };

  auto unite = [&](int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (ufRank[a] < ufRank[b]) std::swap(a, b);
    parent[b] = a;
    if (ufRank[a] == ufRank[b]) ufRank[a]++;
  };

  // A run encodes a maximal consecutive sequence of occupied voxels
  // in the k-direction for a given (i,j).
  struct Run {
    int kStart;   // inclusive
    int kEnd;     // inclusive
    int repIdx;   // global voxel index of the first voxel in the run
  };

  // We keep two slices worth of runs for inter-slice merging.
  // prevSliceRuns[j] = runs from slice (i-1) at row j.
  // currSliceRuns[j] = runs from slice i at row j.
  std::vector<std::vector<Run>> prevSliceRuns(ny);
  std::vector<std::vector<Run>> currSliceRuns(ny);

  // Merge two sorted run-lists from adjacent rows/slices.
  // Two runs overlap in k if their k-ranges share at least one position.
  auto mergeRunLists = [&](const std::vector<Run>& A,
                            const std::vector<Run>& B) {
    std::size_t ai = 0, bi = 0;
    while (ai < A.size() && bi < B.size()) {
      const Run& ra = A[ai];
      const Run& rb = B[bi];
      if (ra.kEnd < rb.kStart) { ++ai; continue; }
      if (rb.kEnd < ra.kStart) { ++bi; continue; }
      // Overlap detected
      unite(ra.repIdx, rb.repIdx);
      // Advance the run that ends earlier; the other may still overlap
      if (ra.kEnd <= rb.kEnd) ++ai;
      else                    ++bi;
    }
  };

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      currSliceRuns[j].clear();

      // Scan row (i,j) along k: collect runs and union voxels within each run
      int kStart = -1;
      int repIdx = -1;

      for (int k = 0; k <= nz; ++k) {
        bool occ = (k < nz) && (occupancy[idx3(i, j, k, ny, nz)] == 1);
        if (occ) {
          int vIdx = static_cast<int>(idx3(i, j, k, ny, nz));
          if (kStart < 0) {
            // Start a new run
            kStart = k;
            repIdx = vIdx;
          } else {
            // Extend existing run: connect this voxel to the run's representative
            unite(repIdx, vIdx);
          }
        } else {
          if (kStart >= 0) {
            // Close the run
            currSliceRuns[j].push_back({kStart, k - 1, find(repIdx)});
            kStart = -1;
            repIdx = -1;
          }
        }
      }

      // Merge with previous row in same slice (j-1 same i)
      if (j > 0 && !currSliceRuns[j - 1].empty()) {
        mergeRunLists(currSliceRuns[j - 1], currSliceRuns[j]);
      }

      // Merge with corresponding row in previous slice (i-1, same j)
      if (i > 0 && !prevSliceRuns[j].empty()) {
        mergeRunLists(prevSliceRuns[j], currSliceRuns[j]);
      }
    }

    // Move current slice runs to previous for next slice iteration
    std::swap(prevSliceRuns, currSliceRuns);
  }

  // Final pass: count clusters
  std::vector<int> clusterSizeMap(totalSize, 0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          int root = find(static_cast<int>(idx));
          clusterSizeMap[root]++;
          if (labels) (*labels)[idx] = root;
          visited[idx] = 1;
          result.visitedVoxels++;
        }
      }
    }
  }

  for (std::size_t i = 0; i < totalSize; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}


// ── Octree-based CCL ─────────────────────────────────────────────────────────
//
// Hierarchical 3D grid subdivision.
// Empty octants are detected early and skipped entirely (O(1) per empty octant).
// At leaf level, local 6-connectivity is resolved via union-find.
// After each recursive split, cross-boundary connections at the three midplanes
// are added by scanning the shared face rows.
//
// Shares BLS's "skip empty space" philosophy, but uses a different strategy:
// crystallographic lattice (BLS) vs. regular octree subdivision (here).

namespace {

// Forward declaration of octree recursive helper
static void octreeRecurse(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<int>& parent,
    std::vector<int>& ufRank,
    int leafSize,
    int x0, int y0, int z0,
    int x1, int y1, int z1);

static int octreeFind(std::vector<int>& parent, int x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

static void octreeUnite(std::vector<int>& parent, std::vector<int>& ufRank,
                         int a, int b) {
  a = octreeFind(parent, a);
  b = octreeFind(parent, b);
  if (a == b) return;
  if (ufRank[a] < ufRank[b]) std::swap(a, b);
  parent[b] = a;
  if (ufRank[a] == ufRank[b]) ufRank[a]++;
}

static void octreeRecurse(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<int>& parent,
    std::vector<int>& ufRank,
    int leafSize,
    int x0, int y0, int z0,
    int x1, int y1, int z1) {

  int dx = x1 - x0;
  int dy = y1 - y0;
  int dz = z1 - z0;
  if (dx <= 0 || dy <= 0 || dz <= 0) return;

  // Early termination: check if octant contains any occupied voxel
  bool hasOccupied = false;
  for (int i = x0; i < x1 && !hasOccupied; ++i)
    for (int j = y0; j < y1 && !hasOccupied; ++j)
      for (int k = z0; k < z1 && !hasOccupied; ++k)
        if (occupancy[idx3(i, j, k, ny, nz)] == 1) hasOccupied = true;

  if (!hasOccupied) return;

  // Leaf: connect all 6-neighbors within this sub-cube
  if (dx <= leafSize && dy <= leafSize && dz <= leafSize) {
    for (int i = x0; i < x1; ++i)
      for (int j = y0; j < y1; ++j)
        for (int k = z0; k < z1; ++k) {
          std::size_t idx = idx3(i, j, k, ny, nz);
          if (occupancy[idx] != 1) continue;
          for (int d = 0; d < 6; ++d) {
            int ni = i + deltas6[d][0];
            int nj = j + deltas6[d][1];
            int nk = k + deltas6[d][2];
            // Only connect within this leaf's bounds
            if (ni >= x0 && ni < x1 && nj >= y0 && nj < y1 && nk >= z0 && nk < z1) {
              std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
              if (occupancy[nIdx] == 1)
                octreeUnite(parent, ufRank, static_cast<int>(idx),
                             static_cast<int>(nIdx));
            }
          }
        }
    return;
  }

  // Split into up to 8 octants at the midpoint of each dimension
  int mx = (dx > 1) ? (x0 + x1) / 2 : x1;
  int my = (dy > 1) ? (y0 + y1) / 2 : y1;
  int mz = (dz > 1) ? (z0 + z1) / 2 : z1;

  // Recurse into 8 children (some may be empty — handled by hasOccupied check above)
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, x0, y0, z0, mx, my, mz);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, x0, y0, mz, mx, my, z1);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, x0, my, z0, mx, y1, mz);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, x0, my, mz, mx, y1, z1);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, mx, y0, z0, x1, my, mz);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, mx, y0, mz, x1, my, z1);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, mx, my, z0, x1, y1, mz);
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize, mx, my, mz, x1, y1, z1);

  // Merge across X midplane (mx-1 | mx)
  if (mx > x0 && mx < x1) {
    for (int j = y0; j < y1; ++j)
      for (int k = z0; k < z1; ++k) {
        std::size_t i1 = idx3(mx - 1, j, k, ny, nz);
        std::size_t i2 = idx3(mx,     j, k, ny, nz);
        if (occupancy[i1] == 1 && occupancy[i2] == 1)
          octreeUnite(parent, ufRank, static_cast<int>(i1), static_cast<int>(i2));
      }
  }

  // Merge across Y midplane (my-1 | my)
  if (my > y0 && my < y1) {
    for (int i = x0; i < x1; ++i)
      for (int k = z0; k < z1; ++k) {
        std::size_t i1 = idx3(i, my - 1, k, ny, nz);
        std::size_t i2 = idx3(i, my,     k, ny, nz);
        if (occupancy[i1] == 1 && occupancy[i2] == 1)
          octreeUnite(parent, ufRank, static_cast<int>(i1), static_cast<int>(i2));
      }
  }

  // Merge across Z midplane (mz-1 | mz)
  if (mz > z0 && mz < z1) {
    for (int i = x0; i < x1; ++i)
      for (int j = y0; j < y1; ++j) {
        std::size_t i1 = idx3(i, j, mz - 1, ny, nz);
        std::size_t i2 = idx3(i, j, mz,     ny, nz);
        if (occupancy[i1] == 1 && occupancy[i2] == 1)
          octreeUnite(parent, ufRank, static_cast<int>(i1), static_cast<int>(i2));
      }
  }
}

}  // anonymous namespace (extended)

ClusterResult octreeCCL(
    int nx, int ny, int nz,
    int leafSize,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  if (leafSize <= 0) leafSize = 8;

  std::vector<int> parent(totalSize);
  std::vector<int> ufRank(totalSize, 0);
  for (std::size_t i = 0; i < totalSize; ++i) parent[i] = static_cast<int>(i);

  // Recursive octree processing
  octreeRecurse(nx, ny, nz, occupancy, parent, ufRank, leafSize,
                0, 0, 0, nx, ny, nz);

  // Count clusters
  std::vector<int> clusterSizeMap(totalSize, 0);
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
        std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] == 1) {
          int root = octreeFind(parent, static_cast<int>(idx));
          clusterSizeMap[root]++;
          if (labels) (*labels)[idx] = root;
          visited[idx] = 1;
          result.visitedVoxels++;
        }
      }

  for (std::size_t i = 0; i < totalSize; ++i) {
    if (clusterSizeMap[i] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizeMap[i]);
      result.maxCluster = std::max(result.maxCluster, clusterSizeMap[i]);
    }
  }

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}


// ── VCCS — Voxel Cloud Connected Segmentation ─────────────────────────────────
//
// Places seed voxels on a uniform 3D lattice (spacing = seedResolution voxels).
// Expands regions from seeds via a Dijkstra-like BFS ordered by Euclidean
// distance to the seed origin. The result is a Voronoi partition of occupied
// voxels anchored at the seed positions.
//
// Key distinction from BLS:
//   - VCCS uses uniform grid seeding → many seeds land in empty space (wasted).
//   - BLS uses crystallographic lattice seeding aligned to the physical structure.
//   - VCCS does NOT guarantee topological correctness: two disconnected but
//     spatially close regions can be assigned to the same supervoxel.
//
// Parameter: seedResolution (--algo-eps, default 3.0 voxels) = seed spacing.

ClusterResult vccs(
    int nx, int ny, int nz,
    double seedResolution,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  int S = std::max(1, static_cast<int>(std::round(seedResolution)));

  // assignment[idx] = cluster id (or -1 = unassigned)
  std::vector<int> assignment(totalSize, -1);

  // Priority queue: (dist_to_seed, voxel_idx, cluster_id, seed_i, seed_j, seed_k)
  // Ordered by ascending distance so nearest-seed claims each voxel first.
  struct PQEntry {
    double dist;
    int idx;
    int clusterId;
    int si, sj, sk;
    bool operator>(const PQEntry& o) const { return dist > o.dist; }
  };
  std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

  int clusterId = 0;

  // Place seeds on a regular 3D grid.
  // Offset by S/2 so seeds are centred within their cells.
  int halfS = S / 2;
  for (int si = halfS; si < nx; si += S) {
    for (int sj = halfS; sj < ny; sj += S) {
      for (int sk = halfS; sk < nz; sk += S) {
        std::size_t seedIdx = idx3(si, sj, sk, ny, nz);
        if (occupancy[seedIdx] == 1 && assignment[seedIdx] == -1) {
          // Start a new cluster from this seed
          assignment[seedIdx] = clusterId;
          visited[seedIdx] = 1;
          result.visitedVoxels++;

          // Push seed's 6-neighbours into the priority queue
          for (int d = 0; d < 6; ++d) {
            int ni = si + deltas6[d][0];
            int nj = sj + deltas6[d][1];
            int nk = sk + deltas6[d][2];
            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
              std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
              if (occupancy[nIdx] == 1 && assignment[nIdx] == -1) {
                double dd = std::sqrt(static_cast<double>((ni - si) * (ni - si) +
                                                           (nj - sj) * (nj - sj) +
                                                           (nk - sk) * (nk - sk)));
                pq.push({dd, static_cast<int>(nIdx), clusterId, si, sj, sk});
              }
            }
          }
          ++clusterId;
        }
      }
    }
  }

  // BFS expansion: assign each unassigned occupied voxel to the nearest seed's cluster
  while (!pq.empty()) {
    PQEntry e = pq.top();
    pq.pop();

    if (assignment[e.idx] != -1) continue;  // already claimed
    if (occupancy[e.idx] != 1) continue;     // sanity check

    assignment[e.idx] = e.clusterId;
    visited[e.idx] = 1;
    result.visitedVoxels++;

    // Decompose flat index back to (i,j,k)
    int nynz = ny * nz;
    int i0 = e.idx / nynz;
    int rem = e.idx % nynz;
    int j0 = rem / nz;
    int k0 = rem % nz;

    for (int d = 0; d < 6; ++d) {
      int ni = i0 + deltas6[d][0];
      int nj = j0 + deltas6[d][1];
      int nk = k0 + deltas6[d][2];
      if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
        std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
        if (occupancy[nIdx] == 1 && assignment[nIdx] == -1) {
          double dd = std::sqrt(static_cast<double>((ni - e.si) * (ni - e.si) +
                                                     (nj - e.sj) * (nj - e.sj) +
                                                     (nk - e.sk) * (nk - e.sk)));
          pq.push({dd, static_cast<int>(nIdx), e.clusterId, e.si, e.sj, e.sk});
        }
      }
    }
  }

  // Any occupied voxels not reached by any seed (seed grid misses them entirely)
  // are assigned as new isolated clusters via simple DFS.
  // This illustrates VCCS's coverage gap vs. BLS's lattice-guaranteed coverage.
  std::vector<StackElement> stack;
  stack.reserve(1000);
  for (std::size_t flatIdx = 0; flatIdx < totalSize; ++flatIdx) {
    if (occupancy[flatIdx] != 1 || assignment[flatIdx] != -1) continue;

    // Start a new cluster from this unreached voxel
    int i0 = static_cast<int>(flatIdx) / (ny * nz);
    int rem = static_cast<int>(flatIdx) % (ny * nz);
    int j0 = rem / nz;
    int k0 = rem % nz;

    stack.clear();
    stack.push_back({i0, j0, k0, 0});
    assignment[flatIdx] = clusterId;

    while (!stack.empty()) {
      StackElement curr = stack.back();
      stack.pop_back();
      std::size_t idx = idx3(curr.i, curr.j, curr.k, ny, nz);
      if (assignment[idx] == clusterId) {
        visited[idx] = 1;
        result.visitedVoxels++;
        for (int d = 0; d < 6; ++d) {
          int ni = curr.i + deltas6[d][0];
          int nj = curr.j + deltas6[d][1];
          int nk = curr.k + deltas6[d][2];
          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            std::size_t nIdx = idx3(ni, nj, nk, ny, nz);
            if (occupancy[nIdx] == 1 && assignment[nIdx] == -1) {
              assignment[nIdx] = clusterId;
              stack.push_back({ni, nj, nk, 0});
            }
          }
        }
      }
    }
    ++clusterId;
  }

  // Tally cluster sizes
  if (clusterId == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  std::vector<int> clusterSizes(clusterId, 0);
  for (std::size_t idx = 0; idx < totalSize; ++idx) {
    if (occupancy[idx] == 1 && assignment[idx] >= 0) {
      clusterSizes[assignment[idx]]++;
    }
  }

  for (int c = 0; c < clusterId; ++c) {
    if (clusterSizes[c] > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(clusterSizes[c]);
      result.maxCluster = std::max(result.maxCluster, clusterSizes[c]);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace bls
