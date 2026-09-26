#include "cluster/Algorithms.hpp"

#include <algorithm>
#include <set>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
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
  if (lower == "gcbd" || lower == "union_find" || lower == "uf") return ClusterAlgorithm::GCBD;
  if (lower == "hdbscan") return ClusterAlgorithm::HDBSCAN;
  if (lower == "cc3d" || lower == "cc3d_basic") return ClusterAlgorithm::CC3D;
  if (lower == "cc3d_optimized" || lower == "cc3d_opt") return ClusterAlgorithm::CC3DOptimized;
  if (lower == "rle_ccl" || lower == "rleccl" || lower == "rle") return ClusterAlgorithm::RLECCL;
  if (lower == "rle_ccl_optimized" || lower == "rleccloptimized")
    return ClusterAlgorithm::RLECCLOptimized;
  if (lower == "vccs") return ClusterAlgorithm::VCCS;
  if (lower == "vccs_optimized" || lower == "vccsoptimized")
    return ClusterAlgorithm::VCCSOptimized;

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
    case ClusterAlgorithm::GCBD: return "gcbd";
    case ClusterAlgorithm::HDBSCAN: return "hdbscan";
    case ClusterAlgorithm::CC3D: return "cc3d";
    case ClusterAlgorithm::CC3DOptimized: return "cc3d_optimized";
    case ClusterAlgorithm::RLECCL: return "rle_ccl";
    case ClusterAlgorithm::RLECCLOptimized: return "rle_ccl_optimized";
    case ClusterAlgorithm::VCCS: return "vccs";
    case ClusterAlgorithm::VCCSOptimized: return "vccs_optimized";
  }
  return "unknown";
}

std::vector<std::string> listAlgorithms() {
  return {"bls", "traditional_dfs", "skip_dfs", "dbscan",
          "hierarchical", "kmeans", "gcbd", "hdbscan",
          "cc3d", "cc3d_optimized", "rle_ccl", "rle_ccl_optimized",
          "vccs", "vccs_optimized"};
}

bool supportsLabels(ClusterAlgorithm algo) {
  switch (algo) {
    case ClusterAlgorithm::TraditionalDFS:
    case ClusterAlgorithm::GCBD:
    case ClusterAlgorithm::CC3D:
    case ClusterAlgorithm::CC3DOptimized:
    case ClusterAlgorithm::RLECCL:
    case ClusterAlgorithm::RLECCLOptimized:
    case ClusterAlgorithm::SkipDFS:       // audit plan item 7: partitions, not counts, in tests
    case ClusterAlgorithm::Hierarchical:
      return true;
    default:
      // BLS is labelled through Analyzer, not here. The remaining methods
      // (dbscan, hdbscan -- withdrawn --, kmeans and both vccs variants) have no
      // label output.
      return false;
  }
}

bool supportsPeriodic(ClusterAlgorithm algo) {
  switch (algo) {
    case ClusterAlgorithm::BLS:  // Analyzer::labelGrid (commensurate lattice, periodic refinement)
    case ClusterAlgorithm::TraditionalDFS:
    case ClusterAlgorithm::SkipDFS:
    case ClusterAlgorithm::GCBD:
    case ClusterAlgorithm::CC3DOptimized:
    case ClusterAlgorithm::RLECCLOptimized:
    case ClusterAlgorithm::VCCSOptimized:
      return true;
    default:
      // Withdrawn textbook tracks (cc3d, rle_ccl, vccs), withdrawn DBSCAN/HDBSCAN (D4/D5),
      // k-means (declared non-periodic, D11), hierarchical (SI only).
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

  if (params.periodic && !supportsPeriodic(algo)) {
    throw std::runtime_error("PBC is not implemented for " + algorithmToString(algo));
  }

  // Note: BLS is handled separately in the main analyzer since it requires
  // lattice enumeration. This function handles the comparison algorithms.
  switch (algo) {
    case ClusterAlgorithm::BLS:
      throw std::runtime_error("BLS algorithm should be run through the Analyzer class.");
    case ClusterAlgorithm::TraditionalDFS:
      return traditionalDFS(params.nx, params.ny, params.nz, occupancy, visited, labels,
                            params.periodic);
    case ClusterAlgorithm::SkipDFS:
      return skipDFS(params.nx, params.ny, params.nz, params.skipDfsJumpDistance, occupancy,
                     visited, labels, params.periodic);
    case ClusterAlgorithm::DBSCAN:
      return dbscan(params.nx, params.ny, params.nz, params.eps, params.minPts, occupancy, visited);
    case ClusterAlgorithm::Hierarchical:
      return hierarchical(params.nx, params.ny, params.nz, params.threshold, occupancy, visited,
                          labels);
    case ClusterAlgorithm::KMeans:
      return kmeans(params.nx, params.ny, params.nz, params.k, occupancy, visited);
    case ClusterAlgorithm::GCBD:
      return gcbd(params.nx, params.ny, params.nz, occupancy, visited, labels, params.periodic);
    case ClusterAlgorithm::HDBSCAN:
      return hdbscan(params.nx, params.ny, params.nz, params.minClusterSize, params.minSamples, occupancy, visited);
    case ClusterAlgorithm::CC3D:
      return cc3d(params.nx, params.ny, params.nz, params.connectivity, occupancy, visited, labels);
    case ClusterAlgorithm::CC3DOptimized:
      return cc3dOptimized(params.nx, params.ny, params.nz, params.connectivity, occupancy, visited,
                           labels, params.periodic);
    case ClusterAlgorithm::RLECCL:
      return rleCCL(params.nx, params.ny, params.nz, occupancy, visited, labels);
    case ClusterAlgorithm::RLECCLOptimized:
      return rleCCLOptimized(params.nx, params.ny, params.nz, occupancy, visited, labels,
                             params.periodic);
    case ClusterAlgorithm::VCCS:
      // Use params.eps as seed spacing in voxels (default 3.0)
      return vccs(params.nx, params.ny, params.nz, params.eps, occupancy, visited);
    case ClusterAlgorithm::VCCSOptimized:
      return vccsOptimized(params.nx, params.ny, params.nz, params.eps, occupancy, visited,
                           params.connectivity, params.periodic);
  }
  throw std::runtime_error("Unhandled algorithm type");
}

// Traditional DFS - kept simple for benchmarking
namespace {

template <bool Periodic>
ClusterResult traditionalDFSImpl(
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
          if constexpr (Periodic) {
            // PBC xyz: the neighbour across a face is the voxel on the opposite face.
            ni = ni < 0 ? ni + nx : (ni >= nx ? ni - nx : ni);
            nj = nj < 0 ? nj + ny : (nj >= ny ? nj - ny : nj);
            nk = nk < 0 ? nk + nz : (nk >= nz ? nk - nz : nk);
            stack.push_back({ni, nj, nk, 0});
          } else {
            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
              stack.push_back({ni, nj, nk, 0});
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

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult traditionalDFS(int nx, int ny, int nz, const std::vector<uint8_t>& occupancy,
                             std::vector<uint8_t>& visited, std::vector<int>* labels, bool periodic) {
  // PBC xyz runs a separate instantiation: PBC none executes the pre-audit code
  // unchanged (same instructions, same cost); only the periodic run pays for wrapping.
  return periodic ? traditionalDFSImpl<true>(nx, ny, nz, occupancy, visited, labels)
                  : traditionalDFSImpl<false>(nx, ny, nz, occupancy, visited, labels);
}

// Skip-DFS - kept simple for benchmarking
namespace {

template <bool Periodic>
ClusterResult skipDFSImpl(
    int nx, int ny, int nz, int skip,
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

  // In-grid test for a neighbour; under PBC xyz the neighbour is wrapped onto the grid
  // instead (true modulo: a jump can exceed a small dimension).
  auto inside = [nx, ny, nz](int& a, int& b, int& c) {
    if constexpr (Periodic) {
      a = ((a % nx) + nx) % nx;
      b = ((b % ny) + ny) % ny;
      c = ((c % nz) + nz) % nz;
      return true;
    } else {
      return a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz;
    }
  };

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
        if (labels) (*labels)[index] = result.nclusters;  // dense: ordinal of this component
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
              if (inside(ii, jj, kk)) {
                stack.push_back({ii, jj, kk, 0});
              }
            }
            // Enqueue the skip-target
            int ni = i + di;
            int nj = j + dj;
            int nk = k + dk;
            if (inside(ni, nj, nk)) {
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
            if (inside(ni, nj, nk)) {
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

  compactLabels(labels);
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult skipDFS(int nx, int ny, int nz, int skip, const std::vector<uint8_t>& occupancy,
                      std::vector<uint8_t>& visited, std::vector<int>* labels, bool periodic) {
  // PBC xyz runs a separate instantiation: PBC none executes the pre-audit code
  // unchanged (same instructions, same cost); only the periodic run pays for wrapping.
  return periodic ? skipDFSImpl<true>(nx, ny, nz, skip, occupancy, visited, labels)
                  : skipDFSImpl<false>(nx, ny, nz, skip, occupancy, visited, labels);
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
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, static_cast<std::size_t>(nx) * ny * nz);

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
    if (labels) (*labels)[occIdx] = root;
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

  compactLabels(labels);
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

  if (k < 1) {
    // KMEANS-1: k <= 0 used to size the centroid arrays to zero and then write into them.
    throw std::invalid_argument("kmeans: k must be >= 1 (got " + std::to_string(k) + ")");
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

  // Lloyd iterations, at most 10 (the campaign's cap, disclosed). Stops early at a fixed
  // point: when no assignment changes, the centroids the next update would compute are
  // the ones already in place, so every later iteration repeats this one -- the result is
  // identical to running all 10, only the redundant work is skipped (KMEANS-2 hygiene, D11).
  for (int iter = 0; iter < 10; ++iter) {
    bool changed = (iter == 0);
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
      if (assignments[p] != bestK) changed = true;
      assignments[p] = bestK;
    }
    if (!changed) break;

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

// GCBD (Grid-based Connectivity using Union-Find)
namespace {

template <bool Periodic>
ClusterResult gcbdImpl(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // In-grid test for a neighbour; under PBC xyz the neighbour is wrapped onto the grid
  // instead (true modulo: a jump can exceed a small dimension).
  auto inside = [nx, ny, nz](int& a, int& b, int& c) {
    if constexpr (Periodic) {
      a = ((a % nx) + nx) % nx;
      b = ((b % ny) + ny) % ny;
      c = ((c % nz) + nz) % nz;
      return true;
    } else {
      return a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz;
    }
  };

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
          if (inside(ni, nj, nk)) {
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

}  // namespace

ClusterResult gcbd(int nx, int ny, int nz, const std::vector<uint8_t>& occupancy,
                   std::vector<uint8_t>& visited, std::vector<int>* labels, bool periodic) {
  // PBC xyz runs a separate instantiation: PBC none executes the pre-audit code
  // unchanged (same instructions, same cost); only the periodic run pays for wrapping.
  return periodic ? gcbdImpl<true>(nx, ny, nz, occupancy, visited, labels)
                  : gcbdImpl<false>(nx, ny, nz, occupancy, visited, labels);
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

// ── CC3D, optimized track ────────────────────────────────────────────────────
//
// The method as it is actually deployed, rather than the textbook form the
// fair track holds to. Two changes, both of which matter here:
//
// 1. SAUF (Scan plus Array-based Union-Find, Wu/Otoo/Suzuki 2005): a two-pass
//    raster scan in which each voxel consults only its ALREADY-SCANNED
//    neighbours -- three of the six in 6-connectivity -- and takes an existing
//    provisional label instead of creating and then merging one. The fair
//    track instead probes all six neighbours of every occupied voxel and
//    unions in both directions, doing every merge twice.
//
// 2. Every auxiliary array is sized to the number of OCCUPIED voxels, not to
//    the grid volume. This is the change that dominates on E1-like data.
//    Task 11 measured the fair cc3d at 81 ms on a completely EMPTY 23.72M-voxel
//    grid -- no occupied voxels, so no clustering work whatsoever -- because it
//    allocates and initialises `parent` and `clusterSizeMap` at grid size. The
//    previous "optimized" variant added a third grid-sized array (`rank`) and
//    was consequently 21-24% SLOWER than the fair one at -O3, which is what
//    manuscript Table 2 reports. Union-find asymptotics were never the problem.
//
// No array here is sized to the grid volume. The backward stencil reaches at
// most one plane back in i, so the scan carries a rolling TWO-PLANE buffer of
// provisional labels (2*ny*nz ints, 660 KB at E1 size against 95 MB for a full
// grid) and records each occupied voxel's provisional label into a compact list
// as it goes, so the resolve pass never needs to address by coordinate again.
namespace {

template <bool Periodic>
ClusterResult cc3dOptimizedImpl(
    int nx, int ny, int nz,
    int connectivity,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // Rolling provisional-label planes: `prevPlane` is i-1, `currPlane` is i.
  // 0 means "no label". This is the whole of the coordinate-addressed state.
  const std::size_t planeYZ = static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
  std::vector<int> prevPlane(planeYZ, 0);
  std::vector<int> currPlane(planeYZ, 0);

  // Union-find over PROVISIONAL LABELS, not voxels, so its size tracks
  // occupancy rather than volume. Index 0 is reserved as "no label".
  std::vector<int> parent;
  parent.reserve(1024);
  parent.push_back(0);

  auto find = [&parent](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];  // path halving
      x = parent[x];
    }
    return x;
  };
  auto unite = [&](int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (a < b) parent[b] = a; else parent[a] = b;  // keep the smaller root
  };

  // Occupied voxels and their provisional labels, in raster order. Sized to
  // occupancy; this is what the resolve pass walks instead of the grid.
  std::vector<int> occIdx, occLab;

  // Backward neighbour offsets: those already visited under an (i,j,k) raster
  // order. Three for 6-connectivity, thirteen for 26 -- exactly half of each
  // stencil, which is the whole point of the scan order. None reaches further
  // back than i-1, which is what makes the two-plane buffer sufficient.
  constexpr int back6[3][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
  constexpr int back26[13][3] = {
      {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
      {-1, -1, 0}, {-1, 1, 0}, {-1, 0, -1}, {-1, 0, 1},
      {0, -1, -1}, {0, -1, 1},
      {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}, {-1, 1, 1}};
  const int (*back)[3] = (connectivity == 26) ? back26 : back6;
  const int numBack   = (connectivity == 26) ? 13 : 3;

  // Pass 1: assign provisional labels, recording equivalences.
  for (int i = 0; i < nx; ++i) {
    std::fill(currPlane.begin(), currPlane.end(), 0);
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t idx = idx3(i, j, k, ny, nz);
        if (occupancy[idx] != 1) continue;

        int lab = 0;
        for (int d = 0; d < numBack; ++d) {
          const int ni = i + back[d][0];
          const int nj = j + back[d][1];
          const int nk = k + back[d][2];
          if (ni < 0 || nj < 0 || nj >= ny || nk < 0 || nk >= nz) continue;
          const std::size_t off = static_cast<std::size_t>(nj) * static_cast<std::size_t>(nz) +
                                  static_cast<std::size_t>(nk);
          const int nl = (ni == i) ? currPlane[off] : prevPlane[off];
          if (nl == 0) continue;
          if (lab == 0) lab = nl;
          else          unite(lab, nl);
        }
        if (lab == 0) {                       // no labelled backward neighbour
          lab = static_cast<int>(parent.size());
          parent.push_back(lab);              // new provisional label
        }
        currPlane[static_cast<std::size_t>(j) * static_cast<std::size_t>(nz) +
                  static_cast<std::size_t>(k)] = lab;
        occIdx.push_back(static_cast<int>(idx));
        occLab.push_back(lab);
      }
    }
    std::swap(prevPlane, currPlane);
  }

  if constexpr (Periodic) {
    // PBC xyz: pass 1 saw every neighbour pair inside the grid; the pairs that wrap
    // across a face are united here, from the voxels on the faces, over the full
    // stencil. A neighbour's provisional label is found in the raster-ordered (so
    // ascending) occupied list.
    auto labelAt = [&](int vidx) {
      const auto it = std::lower_bound(occIdx.begin(), occIdx.end(), vidx);
      return (it != occIdx.end() && *it == vidx)
                 ? occLab[static_cast<std::size_t>(it - occIdx.begin())]
                 : 0;
    };
    const int plane = ny * nz;
    for (std::size_t n = 0; n < occIdx.size(); ++n) {
      const int i = occIdx[n] / plane;
      const int j = (occIdx[n] - i * plane) / nz;
      const int k = occIdx[n] - i * plane - j * nz;
      if (i != 0 && i != nx - 1 && j != 0 && j != ny - 1 && k != 0 && k != nz - 1) continue;
      for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
          for (int dk = -1; dk <= 1; ++dk) {
            const int manhattan = std::abs(di) + std::abs(dj) + std::abs(dk);
            if (manhattan == 0 || (connectivity != 26 && manhattan != 1)) continue;
            int ni = i + di, nj = j + dj, nk = k + dk;
            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) continue;
            ni = (ni + nx) % nx;
            nj = (nj + ny) % ny;
            nk = (nk + nz) % nz;
            const int nidx = static_cast<int>(idx3(ni, nj, nk, ny, nz));
            if (occupancy[static_cast<std::size_t>(nidx)] != 1) continue;
            unite(occLab[n], labelAt(nidx));
          }
        }
      }
    }
  }

  // Pass 2: resolve equivalences to dense ids and tally, walking the compact
  // occupied list. Ids are handed out in order of first appearance, which is
  // raster order -- exactly what compactLabels() would produce.
  std::vector<int> finalId(parent.size(), -1);
  std::vector<int> sizes;
  sizes.reserve(parent.size());

  for (std::size_t n = 0; n < occIdx.size(); ++n) {
    const int root = find(occLab[n]);
    int id = finalId[static_cast<std::size_t>(root)];
    if (id < 0) {
      id = static_cast<int>(sizes.size());
      finalId[static_cast<std::size_t>(root)] = id;
      sizes.push_back(0);
    }
    ++sizes[static_cast<std::size_t>(id)];
    const std::size_t idx = static_cast<std::size_t>(occIdx[n]);
    if (labels) (*labels)[idx] = id;
    visited[idx] = 1;
    result.visitedVoxels++;
  }

  result.nclusters = static_cast<int>(sizes.size());
  for (int s : sizes) {
    result.clusterSizes.push_back(s);
    result.maxCluster = std::max(result.maxCluster, s);
  }

  // Labels are already dense and assigned in order of first appearance, which
  // is exactly what compactLabels() would produce, so it is not called here.
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult cc3dOptimized(int nx, int ny, int nz, int connectivity,
                            const std::vector<uint8_t>& occupancy, std::vector<uint8_t>& visited,
                            std::vector<int>* labels, bool periodic) {
  // PBC xyz runs a separate instantiation: PBC none executes the pre-audit code
  // unchanged (same instructions, same cost); only the periodic run pays for wrapping.
  return periodic ? cc3dOptimizedImpl<true>(nx, ny, nz, connectivity, occupancy, visited, labels)
                  : cc3dOptimizedImpl<false>(nx, ny, nz, connectivity, occupancy, visited, labels);
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



// ── RLE-CCL, optimized track ─────────────────────────────────────────────────
//
// He, Chao, Suzuki & Wu, "A Run-Based Two-Scan Labeling Algorithm",
// IEEE Trans. Image Processing 17(5), 2008 -- the paper the fair track already
// cites but does not follow.
//
// Three things separate this from the fair track, all of them the point of RLE:
//
// 1. RUNS, not voxels, are the union-find domain. The fair track builds runs
//    and then unions every voxel inside each run to the run's representative
//    (Algorithms.cpp, rleCCL, "Extend existing run"), which is precisely the
//    work run-length encoding exists to avoid: a run of length L costs L-1
//    unions there and zero here. The union-find array is sized to the number of
//    runs, which on E1-like data is a few thousand against 23.7M voxels.
//
// 2. Rows are iterated sparsely. The fair track's inner loop visits every k in
//    [0,nz) for every (i,j) -- a full O(NX*NY*NZ) sweep -- and then sweeps the
//    grid twice more to count. Here the row scan is the only pass over the
//    occupancy array, and the counting pass walks the run table.
//
// 3. Adjacency between runs is resolved by the two-pointer merge of sorted run
//    lists (as in the fair track) but applied to run ids, so a merge costs one
//    union per overlapping PAIR rather than one per shared voxel.
//
// Output is identical to the fair track by construction: both compute the
// 6-connected components of the same occupancy, and connectivity of a run to
// its neighbours is unchanged by whether its interior was unioned voxel by
// voxel.
namespace {

template <bool Periodic>
ClusterResult rleCCLOptimizedImpl(
    int nx, int ny, int nz,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited,
    std::vector<int>* labels) {

  ScopedTimer timer;
  ClusterResult result;

  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);
  initLabels(labels, totalSize);

  // Run table. Sized to the number of runs, which is data-dependent.
  struct Run { int i, j, kStart, kEnd; };
  std::vector<Run> runs;
  runs.reserve(1024);

  // Union-find over RUN IDS.
  std::vector<int> parent;
  parent.reserve(1024);
  auto find = [&parent](int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
  };
  auto unite = [&](int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (a < b) parent[b] = a; else parent[a] = b;
  };

  // Run ids for the current and previous row within a slice, and for the whole
  // previous slice. Runs in each list are ordered by kStart, which the scan
  // produces for free and the two-pointer merge below relies on.
  std::vector<std::vector<int>> prevSlice(ny), currSlice(ny);

  // Two-pointer merge of two kStart-sorted run lists: one union per overlapping
  // pair. Runs overlap when their k-ranges share at least one position.
  auto mergeRuns = [&](const std::vector<int>& A, const std::vector<int>& B) {
    std::size_t ai = 0, bi = 0;
    while (ai < A.size() && bi < B.size()) {
      const Run& ra = runs[static_cast<std::size_t>(A[ai])];
      const Run& rb = runs[static_cast<std::size_t>(B[bi])];
      if (ra.kEnd < rb.kStart) { ++ai; continue; }
      if (rb.kEnd < ra.kStart) { ++bi; continue; }
      unite(A[ai], B[bi]);
      if (ra.kEnd <= rb.kEnd) ++ai; else ++bi;
    }
  };

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      currSlice[j].clear();

      // The single pass over the occupancy array. Runs are emitted whole; no
      // per-voxel union-find operation happens anywhere in this function.
      const std::size_t rowBase = idx3(i, j, 0, ny, nz);
      int kStart = -1;
      for (int k = 0; k <= nz; ++k) {
        const bool occ = (k < nz) && (occupancy[rowBase + static_cast<std::size_t>(k)] == 1);
        if (occ) {
          if (kStart < 0) kStart = k;
        } else if (kStart >= 0) {
          const int id = static_cast<int>(runs.size());
          runs.push_back(Run{i, j, kStart, k - 1});
          parent.push_back(id);
          currSlice[j].push_back(id);
          kStart = -1;
        }
      }

      if (j > 0 && !currSlice[j - 1].empty()) mergeRuns(currSlice[j - 1], currSlice[j]);
      if (i > 0 && !prevSlice[j].empty())     mergeRuns(prevSlice[j],     currSlice[j]);
    }
    std::swap(prevSlice, currSlice);
  }

  if constexpr (Periodic) {
    // PBC xyz: runs that touch opposite faces are merged here, after the in-grid
    // merges above. Runs were created in raster order, so each row's runs form a
    // contiguous, kStart-ordered id range.
    std::vector<int> rowStart(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) + 1, 0);
    for (const Run& r : runs)
      ++rowStart[static_cast<std::size_t>(r.i) * static_cast<std::size_t>(ny) +
                 static_cast<std::size_t>(r.j) + 1];
    for (std::size_t r = 1; r < rowStart.size(); ++r) rowStart[r] += rowStart[r - 1];
    auto rowRuns = [&](int i, int j) {
      const std::size_t row = static_cast<std::size_t>(i) * static_cast<std::size_t>(ny) +
                              static_cast<std::size_t>(j);
      std::vector<int> ids;
      for (int id = rowStart[row]; id < rowStart[row + 1]; ++id) ids.push_back(id);
      return ids;
    };
    for (std::size_t row = 0; row + 1 < rowStart.size(); ++row) {  // z faces, within a row
      const int first = rowStart[row], last = rowStart[row + 1] - 1;
      if (first <= last && runs[static_cast<std::size_t>(first)].kStart == 0 &&
          runs[static_cast<std::size_t>(last)].kEnd == nz - 1) {
        unite(first, last);
      }
    }
    for (int i = 0; i < nx; ++i) mergeRuns(rowRuns(i, 0), rowRuns(i, ny - 1));  // y faces
    for (int j = 0; j < ny; ++j) mergeRuns(rowRuns(0, j), rowRuns(nx - 1, j));  // x faces
  }

  // Tally over the RUN TABLE, not the grid. Dense ids are assigned in order of
  // first appearance scanning runs in creation order, which is raster order, so
  // this matches what the fair track's compactLabels() produces.
  std::vector<int> finalId(runs.size(), -1);
  std::vector<int> sizes;
  sizes.reserve(runs.size());

  for (std::size_t r = 0; r < runs.size(); ++r) {
    const int root = find(static_cast<int>(r));
    int id = finalId[static_cast<std::size_t>(root)];
    if (id < 0) {
      id = static_cast<int>(sizes.size());
      finalId[static_cast<std::size_t>(root)] = id;
      sizes.push_back(0);
    }
    const Run& run = runs[r];
    const int len = run.kEnd - run.kStart + 1;
    sizes[static_cast<std::size_t>(id)] += len;
    const std::size_t base = idx3(run.i, run.j, 0, ny, nz);
    for (int k = run.kStart; k <= run.kEnd; ++k) {
      const std::size_t idx = base + static_cast<std::size_t>(k);
      if (labels) (*labels)[idx] = id;
      visited[idx] = 1;
    }
    result.visitedVoxels += static_cast<std::size_t>(len);
  }

  result.nclusters = static_cast<int>(sizes.size());
  for (int sz : sizes) {
    result.clusterSizes.push_back(sz);
    result.maxCluster = std::max(result.maxCluster, sz);
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult rleCCLOptimized(int nx, int ny, int nz, const std::vector<uint8_t>& occupancy,
                              std::vector<uint8_t>& visited, std::vector<int>* labels, bool periodic) {
  // PBC xyz runs a separate instantiation: PBC none executes the pre-audit code
  // unchanged (same instructions, same cost); only the periodic run pays for wrapping.
  return periodic ? rleCCLOptimizedImpl<true>(nx, ny, nz, occupancy, visited, labels)
                  : rleCCLOptimizedImpl<false>(nx, ny, nz, occupancy, visited, labels);
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
        ++result.seedCandidates;
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
  // Captured before the remainder pass below, which also advances clusterId.
  result.seedsPlaced = clusterId;
  result.seedPruneThreshold = 0;   // the fair track does not prune

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


namespace {

// ── VCCS, optimized track ────────────────────────────────────────────────────
//
// PORTED rather than linked: PCL 1.14 pulls in Boost, Eigen, FLANN and VTK for one
// comparison algorithm, this build has no mandatory external dependency, and PCL's entry
// point takes a PointCloud<PointXYZRGBA> with estimated normals -- the adapter from a
// binary voxel grid would be more code than the algorithm. No grid-sized array is
// allocated: the occupied voxels are collected once into a raster-ordered list and every
// per-voxel state is indexed by position in it (voxel index -> position is a binary search).
// PCL's octree adjacency graph is not ported either: on a regular grid a voxel's
// neighbours are index arithmetic away.
//
// VCCS in its published form (audit decision D6; findings VCCS-1, -2, -5, -6):
// a port of PCL's SupervoxelClustering (Papon et al., CVPR 2013; PCL
// segmentation/impl/supervoxel_clustering.hpp, BSD-3-Clause, see THIRD_PARTY_NOTICES.md)
// onto the occupancy grid, with the voxel resolution = 1 voxel and the spatial
// distance term only (every voxel here is identical apart from its position, so
// PCL's colour and normal terms carry no information; manuscript l.446-450).
//
//   seeding   one candidate per occupied seed cell (edge S voxels): the occupied voxel
//             nearest the cell centre, over ALL occupied voxels (PCL: kd-tree nearest to
//             the seed-octree leaf centre); candidates chosen twice are kept once.
//   pruning   keep a candidate iff  #{occupied v : |v - c|^2 < R^2} > 0.05*pi*R^2,
//             R = S/2 (PCL selectInitialSupervoxelSeeds; strict < as FLANN's radius search).
//   growth    max_depth = int(1.8*S); for i = 1 .. max_depth-1: every supervoxel, in seed
//             order, scans the neighbours (PCL: 26; here the method's connectivity, see
//             below) of all its voxels and takes a neighbour whose
//             distance to the supervoxel's CENTROID is below the neighbour's recorded
//             distance (taking it from its previous owner); then empty supervoxels are
//             dropped and centroids recomputed (PCL expandSupervoxels / SupervoxelHelper::
//             expand / updateCentroid). As in PCL, a seed voxel's recorded distance starts
//             at +inf, so a neighbouring supervoxel may take it.
//   output    one cluster per non-empty supervoxel. Voxels no supervoxel reached (PCL leaves
//             them unlabelled) are promoted, one cluster per connected remainder, so that
//             every method is tallied on the same terms (a coverage gap shows up as extra
//             clusters, not as missing voxels).
//
// Deliberate differences from PCL, all forced by determinism or the voxel grid: voxels of a
// supervoxel are visited in voxel-index order (PCL: a pointer-ordered std::set, i.e. memory
// layout); seed cells are aligned with the grid, [cS, cS+S-1] (PCL: the seed octree's
// bounding box, itself an arbitrary alignment); distances are in voxels, not scaled by spatial_importance/S
// (a common positive factor, so every comparison is unchanged).
// Under PBC xyz (Periodic) adjacency wraps, and every distance and centroid is taken with
// the minimum image; the seed-cell grid is not periodic (PCL has no PBC; a cell cut by a face
// is simply a smaller cell).
template <bool Periodic>
ClusterResult vccsOptimizedImpl(int nx, int ny, int nz, double seedResolution, int connectivity,
                                const std::vector<uint8_t>& occupancy,
                                std::vector<uint8_t>& visited) {
  ScopedTimer timer;
  ClusterResult result;

  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  const int S = std::max(1, static_cast<int>(std::round(seedResolution)));
  // Growth / remainder adjacency: PCL's is 26 (face, edge, vertex). Here it is the
  // method's connectivity (6, 18 or 26) so the caller decides whether VCCS refines the
  // same components the exact methods label (6) or PCL's (26); see decision D6.
  std::vector<std::array<int, 3>> nbr;
  for (int di = -1; di <= 1; ++di)
    for (int dj = -1; dj <= 1; ++dj)
      for (int dk = -1; dk <= 1; ++dk) {
        const int m = std::abs(di) + std::abs(dj) + std::abs(dk);
        if (m == 0 || (connectivity == 6 && m != 1) || (connectivity == 18 && m == 3)) continue;
        nbr.push_back({di, dj, dk});
      }

  // The one pass over the occupancy array. Raster order, so `occVox` comes out
  // sorted and can be binary-searched without an explicit sort.
  std::vector<int> occVox;
  for (std::size_t i = 0; i < totalSize; ++i) {
    if (occupancy[i] == 1) occVox.push_back(static_cast<int>(i));
  }
  const int nOcc = static_cast<int>(occVox.size());
  if (nOcc == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }

  auto compactOf = [&occVox](int voxIdx) -> int {
    auto it = std::lower_bound(occVox.begin(), occVox.end(), voxIdx);
    if (it == occVox.end() || *it != voxIdx) return -1;
    return static_cast<int>(it - occVox.begin());
  };
  const int planeYZ = ny * nz;
  auto decode = [&](int voxIdx, int& i, int& j, int& k) {
    i = voxIdx / planeYZ;
    const int rem = voxIdx - i * planeYZ;
    j = rem / nz;
    k = rem % nz;
  };
  // Neighbour / offset voxel -> compact index, -1 when outside (non-periodic) or empty.
  auto occupiedAt = [&](int i, int j, int k) -> int {
    if constexpr (Periodic) {
      i = ((i % nx) + nx) % nx;
      j = ((j % ny) + ny) % ny;
      k = ((k % nz) + nz) % nz;
    } else {
      if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return -1;
    }
    const int v = static_cast<int>(idx3(i, j, k, ny, nz));
    return occupancy[static_cast<std::size_t>(v)] == 1 ? compactOf(v) : -1;
  };
  // Displacement a - b along one axis of length n (minimum image under PBC).
  auto delta = [](double a, double b, int n) {
    double d = a - b;
    if constexpr (Periodic) d -= n * std::round(d / n);
    return d;
  };
  auto dist2 = [&](double ai, double aj, double ak, double bi, double bj, double bk) {
    const double di = delta(ai, bi, nx), dj = delta(aj, bj, ny), dk = delta(ak, bk, nz);
    return di * di + dj * dj + dk * dk;
  };

  // --- 1. seeding: nearest occupied voxel to each occupied seed-cell centre ------------
  // Seed cells are the S^3 blocks [cS, cS+S-1] per axis, centre cS + (S-1)/2 in voxel-
  // centre coordinates. The best voxel inside the cell bounds the search: only voxels
  // strictly closer to the centre (or as close, with a smaller index) can replace it,
  // and they lie in the cube of that radius around the centre.
  std::vector<std::pair<std::array<int, 3>, int>> byCell;
  byCell.reserve(static_cast<std::size_t>(nOcc));
  for (int c = 0; c < nOcc; ++c) {
    int i, j, k; decode(occVox[static_cast<std::size_t>(c)], i, j, k);
    byCell.push_back({{i / S, j / S, k / S}, c});
  }
  std::sort(byCell.begin(), byCell.end());
  std::vector<int> candidates;
  const double half = 0.5 * (S - 1);
  for (std::size_t g = 0; g < byCell.size();) {
    const auto cell = byCell[g].first;
    const double cx = cell[0] * S + half, cy = cell[1] * S + half, cz = cell[2] * S + half;
    int best = -1;
    double bestD2 = 0.0;
    for (; g < byCell.size() && byCell[g].first == cell; ++g) {
      const int c = byCell[g].second;
      int i, j, k; decode(occVox[static_cast<std::size_t>(c)], i, j, k);
      const double d2 = dist2(i, j, k, cx, cy, cz);
      if (best < 0 || d2 < bestD2 || (d2 == bestD2 && c < best)) { best = c; bestD2 = d2; }
    }
    const int reach = static_cast<int>(std::ceil(std::sqrt(bestD2))) + 1;
    const int ci = static_cast<int>(std::floor(cx)), cj = static_cast<int>(std::floor(cy)),
              ck = static_cast<int>(std::floor(cz));
    for (int i = ci - reach; i <= ci + reach + 1; ++i)
      for (int j = cj - reach; j <= cj + reach + 1; ++j)
        for (int k = ck - reach; k <= ck + reach + 1; ++k) {
          const int c = occupiedAt(i, j, k);
          if (c < 0 || c == best) continue;
          int vi, vj, vk; decode(occVox[static_cast<std::size_t>(c)], vi, vj, vk);
          const double d2 = dist2(vi, vj, vk, cx, cy, cz);
          if (d2 < bestD2 || (d2 == bestD2 && c < best)) { best = c; bestD2 = d2; }
        }
    candidates.push_back(best);
  }
  // Two cells can pick the same voxel; PCL would then give one voxel to two helpers
  // (the second addLeaf wins the owner). One supervoxel per voxel is kept.
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

  // --- 2. pruning (PCL selectInitialSupervoxelSeeds) ----------------------------------
  const double rSearch = 0.5 * S;
  const double minPoints = 0.05 * rSearch * rSearch * 3.1415926536;
  const int reachR = static_cast<int>(std::ceil(rSearch));
  std::vector<int> seeds;
  seeds.reserve(candidates.size());
  for (int c : candidates) {
    int i, j, k; decode(occVox[static_cast<std::size_t>(c)], i, j, k);
    int num = 0;
    for (int di = -reachR; di <= reachR; ++di)
      for (int dj = -reachR; dj <= reachR; ++dj)
        for (int dk = -reachR; dk <= reachR; ++dk) {
          if (di * di + dj * dj + dk * dk >= rSearch * rSearch) continue;  // strict, as FLANN
          if (occupiedAt(i + di, j + dj, k + dk) >= 0) ++num;
        }
    if (num > minPoints) seeds.push_back(c);
  }
  result.seedCandidates = static_cast<int>(candidates.size());
  result.seedsPlaced = static_cast<int>(seeds.size());
  result.seedPruneThreshold = minPoints;

  // --- 3. growth (PCL expandSupervoxels) -----------------------------------------------
  const int nSv = static_cast<int>(seeds.size());
  std::vector<int> owner(static_cast<std::size_t>(nOcc), -1);
  std::vector<double> distance(static_cast<std::size_t>(nOcc), std::numeric_limits<double>::max());
  std::vector<std::set<int>> leaves(static_cast<std::size_t>(nSv));
  std::vector<std::array<double, 3>> centroid(static_cast<std::size_t>(nSv));
  std::vector<char> alive(static_cast<std::size_t>(nSv), 1);
  for (int s = 0; s < nSv; ++s) {
    const int c = seeds[static_cast<std::size_t>(s)];
    leaves[static_cast<std::size_t>(s)].insert(c);  // addLeaf: owner set, distance left at +inf
    owner[static_cast<std::size_t>(c)] = s;
    int i, j, k; decode(occVox[static_cast<std::size_t>(c)], i, j, k);
    centroid[static_cast<std::size_t>(s)] = {double(i), double(j), double(k)};
  }

  const int maxDepth = static_cast<int>(1.8 * S);
  std::vector<int> newOwned;
  for (int round = 1; round < maxDepth; ++round) {
    for (int s = 0; s < nSv; ++s) {
      if (!alive[static_cast<std::size_t>(s)]) continue;
      newOwned.clear();
      const auto& cen = centroid[static_cast<std::size_t>(s)];
      for (int leaf : leaves[static_cast<std::size_t>(s)]) {
        int i, j, k; decode(occVox[static_cast<std::size_t>(leaf)], i, j, k);
        for (const auto& o : nbr) {
              const int nc = occupiedAt(i + o[0], j + o[1], k + o[2]);
              if (nc < 0 || owner[static_cast<std::size_t>(nc)] == s) continue;
              int vi, vj, vk; decode(occVox[static_cast<std::size_t>(nc)], vi, vj, vk);
              const double d = std::sqrt(dist2(vi, vj, vk, cen[0], cen[1], cen[2]));
              if (d < distance[static_cast<std::size_t>(nc)]) {
                distance[static_cast<std::size_t>(nc)] = d;
                const int prev = owner[static_cast<std::size_t>(nc)];
                if (prev >= 0) leaves[static_cast<std::size_t>(prev)].erase(nc);
                owner[static_cast<std::size_t>(nc)] = s;
                newOwned.push_back(nc);
              }
            }
      }
      leaves[static_cast<std::size_t>(s)].insert(newOwned.begin(), newOwned.end());
    }
    for (int s = 0; s < nSv; ++s) {
      if (!alive[static_cast<std::size_t>(s)]) continue;
      auto& L = leaves[static_cast<std::size_t>(s)];
      if (L.empty()) { alive[static_cast<std::size_t>(s)] = 0; continue; }
      // Centroid of the supervoxel's voxels, unwrapped around its current centroid.
      const auto c0 = centroid[static_cast<std::size_t>(s)];
      double sx = 0, sy = 0, sz = 0;
      for (int leaf : L) {
        int i, j, k; decode(occVox[static_cast<std::size_t>(leaf)], i, j, k);
        sx += c0[0] + delta(i, c0[0], nx);
        sy += c0[1] + delta(j, c0[1], ny);
        sz += c0[2] + delta(k, c0[2], nz);
      }
      const double n = static_cast<double>(L.size());
      centroid[static_cast<std::size_t>(s)] = {sx / n, sy / n, sz / n};
    }
  }

  // --- 4. tally; unreached voxels promoted per connected remainder ---------------------
  std::vector<int> assignment(static_cast<std::size_t>(nOcc), -1);
  int clusterId = 0;
  for (int s = 0; s < nSv; ++s) {
    if (!alive[static_cast<std::size_t>(s)] || leaves[static_cast<std::size_t>(s)].empty()) continue;
    for (int leaf : leaves[static_cast<std::size_t>(s)]) assignment[static_cast<std::size_t>(leaf)] = clusterId;
    ++clusterId;
  }
  std::vector<int> stack;
  for (int c = 0; c < nOcc; ++c) {
    if (assignment[static_cast<std::size_t>(c)] != -1) continue;
    assignment[static_cast<std::size_t>(c)] = clusterId;
    stack.clear();
    stack.push_back(c);
    while (!stack.empty()) {
      const int cur = stack.back();
      stack.pop_back();
      int i0, j0, k0; decode(occVox[static_cast<std::size_t>(cur)], i0, j0, k0);
      for (const auto& o : nbr) {
            const int nc = occupiedAt(i0 + o[0], j0 + o[1], k0 + o[2]);
            if (nc < 0 || assignment[static_cast<std::size_t>(nc)] != -1) continue;
            assignment[static_cast<std::size_t>(nc)] = clusterId;
            stack.push_back(nc);
          }
    }
    ++clusterId;
  }

  std::vector<int> sizes(static_cast<std::size_t>(clusterId), 0);
  for (int c = 0; c < nOcc; ++c) {
    ++sizes[static_cast<std::size_t>(assignment[static_cast<std::size_t>(c)])];
    visited[static_cast<std::size_t>(occVox[static_cast<std::size_t>(c)])] = 1;
    result.visitedVoxels++;
  }
  for (int sz : sizes) {
    if (sz > 0) {
      result.nclusters++;
      result.clusterSizes.push_back(sz);
      result.maxCluster = std::max(result.maxCluster, sz);
    }
  }

  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult vccsOptimized(int nx, int ny, int nz, double seedResolution,
                            const std::vector<uint8_t>& occupancy, std::vector<uint8_t>& visited,
                            int connectivity, bool periodic) {
  return periodic ? vccsOptimizedImpl<true>(nx, ny, nz, seedResolution, connectivity, occupancy, visited)
                  : vccsOptimizedImpl<false>(nx, ny, nz, seedResolution, connectivity, occupancy, visited);
}

}  // namespace bls
