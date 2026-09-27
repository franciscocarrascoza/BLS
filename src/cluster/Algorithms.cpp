#include "cluster/Algorithms.hpp"

#include <algorithm>
#include <set>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
    case ClusterAlgorithm::DBSCAN:  // Ester et al. 1996 with minimum-image distances
    case ClusterAlgorithm::HDBSCAN: // HDBSCAN* (Campello et al. 2013) with minimum-image distances
      return true;
    default:
      // Withdrawn textbook tracks (cc3d, rle_ccl, vccs),
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
      return dbscan(params.nx, params.ny, params.nz, params.eps, params.minPts, occupancy, visited,
                    params.periodic);
    case ClusterAlgorithm::Hierarchical:
      return hierarchical(params.nx, params.ny, params.nz, params.threshold, occupancy, visited,
                          labels);
    case ClusterAlgorithm::KMeans:
      return kmeans(params.nx, params.ny, params.nz, params.k, occupancy, visited);
    case ClusterAlgorithm::GCBD:
      return gcbd(params.nx, params.ny, params.nz, occupancy, visited, labels, params.periodic);
    case ClusterAlgorithm::HDBSCAN:
      return hdbscan(params.nx, params.ny, params.nz, params.minClusterSize, params.minSamples,
                     occupancy, visited, params.periodic);
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
                           params.periodic);
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

namespace {

// ── DBSCAN ────────────────────────────────────────────────────────────────────
//
// As PUBLISHED: M. Ester, H.-P. Kriegel, J. Sander, X. Xu, "A density-based algorithm for
// discovering clusters in large spatial databases with noise", Proc. KDD-96, pp. 226-231,
// Definitions 1-6 and §4.1 (functions DBSCAN and ExpandCluster). Audit 27.09.26
// (decision D4/D5 revised): replaces an expansion that grew clusters from every point
// reached, border points included (finding S1-10).
//
//   N_Eps(p) = { q : dist(p, q) <= Eps } (Definition 1), Euclidean distance between voxel
//   centres in voxel units, p included. SetOfPoints = the occupied voxels in raster order.
//   DBSCAN: for each point in that order, if it is UNCLASSIFIED, ExpandCluster(point, ClId)
//   and take the next ClId when it returns true.
//   ExpandCluster: seeds = regionQuery(point, Eps); fewer than MinPts -> point is NOISE,
//   return false. Otherwise every seed gets ClId (changeClIds), the point is deleted from
//   the seeds, and while seeds remain: currentP = first seed; result = regionQuery(currentP);
//   if |result| >= MinPts, every result point that is UNCLASSIFIED or NOISE gets ClId and,
//   if it was UNCLASSIFIED, is appended to the seeds; delete currentP.
// Noise points belong to no cluster and are not counted. The region query enumerates the
// integer offsets within Eps (the paper assumes a spatial access method; the neighbourhood
// is the same set); results and seeds follow offset order (dx, dy, dz) lexicographically.
// A border point that lies in two clusters "will be assigned to the cluster discovered first"
// (§4.1 text). The pseudocode's changeClIds(seeds, ClId) would, read literally, move such a
// point to the later cluster when it lies in that cluster's first neighbourhood; the stated
// outcome is followed, i.e. changeClIds relabels only points not yet in a cluster (the only
// difference between the two readings; measured on the audit tests: 5 of 61 grids).
// Under PBC xyz, coordinates wrap and distances are minimum-image.
template <bool Periodic>
ClusterResult dbscanImpl(
    int nx, int ny, int nz,
    double eps, int minPts,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  std::vector<int> occVox;
  for (std::size_t i = 0; i < totalSize; ++i)
    if (occupancy[i] == 1) occVox.push_back(static_cast<int>(i));
  const int nOcc = static_cast<int>(occVox.size());
  if (nOcc == 0) {
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }
  auto compactOf = [&occVox](int v) {
    auto it = std::lower_bound(occVox.begin(), occVox.end(), v);
    return (it != occVox.end() && *it == v) ? static_cast<int>(it - occVox.begin()) : -1;
  };
  const int planeYZ = ny * nz;

  // Integer offsets with dist <= Eps, lexicographic.
  const int r = static_cast<int>(std::floor(eps));
  std::vector<std::array<int, 3>> ball;
  for (int di = -r; di <= r; ++di)
    for (int dj = -r; dj <= r; ++dj)
      for (int dk = -r; dk <= r; ++dk)
        if (static_cast<double>(di * di + dj * dj + dk * dk) <= eps * eps) ball.push_back({di, dj, dk});
  // In a periodic axis shorter than the ball's diameter two offsets reach one voxel.
  const bool dedupe = Periodic && (nx < 2 * r + 1 || ny < 2 * r + 1 || nz < 2 * r + 1);

  auto regionQuery = [&](int c, std::vector<int>& out) {
    out.clear();
    const int v = occVox[static_cast<std::size_t>(c)];
    const int i = v / planeYZ, j = (v - i * planeYZ) / nz, k = v % nz;
    for (const auto& o : ball) {
      int a = i + o[0], b = j + o[1], d = k + o[2];
      if constexpr (Periodic) {
        a = ((a % nx) + nx) % nx;
        b = ((b % ny) + ny) % ny;
        d = ((d % nz) + nz) % nz;
      } else {
        if (a < 0 || a >= nx || b < 0 || b >= ny || d < 0 || d >= nz) continue;
      }
      const std::size_t w = idx3(a, b, d, ny, nz);
      if (occupancy[w] == 1) out.push_back(compactOf(static_cast<int>(w)));
    }
    if (dedupe) {
      std::vector<int> seen;
      std::vector<int> uniq;
      for (int q : out)
        if (std::find(seen.begin(), seen.end(), q) == seen.end()) { seen.push_back(q); uniq.push_back(q); }
      out.swap(uniq);
    }
  };

  constexpr int UNCLASSIFIED = -2, NOISE = -1;
  std::vector<int> clId(static_cast<std::size_t>(nOcc), UNCLASSIFIED);
  std::vector<int> seeds, res;

  auto expandCluster = [&](int point, int id) {
    regionQuery(point, seeds);
    if (static_cast<int>(seeds.size()) < minPts) {
      clId[static_cast<std::size_t>(point)] = NOISE;
      return false;
    }
    for (int s : seeds) {  // changeClIds(seeds, ClId), for points not yet in a cluster (see header)
      int& c = clId[static_cast<std::size_t>(s)];
      if (c == UNCLASSIFIED || c == NOISE) c = id;
    }
    seeds.erase(std::remove(seeds.begin(), seeds.end(), point), seeds.end());
    for (std::size_t f = 0; f < seeds.size(); ++f) {  // seeds.first() ... seeds.delete(currentP)
      regionQuery(seeds[f], res);
      if (static_cast<int>(res.size()) >= minPts) {
        for (int q : res) {
          int& c = clId[static_cast<std::size_t>(q)];
          if (c == UNCLASSIFIED || c == NOISE) {
            if (c == UNCLASSIFIED) seeds.push_back(q);
            c = id;
          }
        }
      }
    }
    return true;
  };

  int clusterId = 0;
  for (int p = 0; p < nOcc; ++p)
    if (clId[static_cast<std::size_t>(p)] == UNCLASSIFIED && expandCluster(p, clusterId)) ++clusterId;

  std::vector<int> sizes(static_cast<std::size_t>(clusterId), 0);
  for (int p = 0; p < nOcc; ++p) {
    const int c = clId[static_cast<std::size_t>(p)];
    if (c < 0) continue;
    ++sizes[static_cast<std::size_t>(c)];
    visited[static_cast<std::size_t>(occVox[static_cast<std::size_t>(p)])] = 1;
    result.visitedVoxels++;
  }
  for (int s : sizes) {
    if (s == 0) continue;
    result.nclusters++;
    result.clusterSizes.push_back(s);
    result.maxCluster = std::max(result.maxCluster, s);
  }
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult dbscan(int nx, int ny, int nz, double eps, int minPts,
                     const std::vector<uint8_t>& occupancy, std::vector<uint8_t>& visited,
                     bool periodic) {
  return periodic ? dbscanImpl<true>(nx, ny, nz, eps, minPts, occupancy, visited)
                  : dbscanImpl<false>(nx, ny, nz, eps, minPts, occupancy, visited);
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

  // Lloyd's iteration (Lloyd, IEEE TIT 28(2):129-137, 1982, Method I: "imposing conditions
  // (14) and (17) alternately" -- centre of mass of each set, nearest-quantum partition --
  // until the sequence converges), in 3-D. It stops at the fixed point: when no assignment
  // changes, the next update would reproduce the current centroids. Audit 27.09.26: the
  // campaign's fixed cap of 10 iterations is not in the source and is removed. Lloyd leaves
  // the empty-set case to "appropriate modifications": an empty cluster keeps its centroid.
  // Initial centroids at evenly spaced points of the raster-ordered list (Lloyd: arbitrary).
  const int kMaxIterations = 1000000;  // guard only; convergence is finite for finite data
  for (int iter = 0;; ++iter) {
    if (iter == kMaxIterations) {
      throw std::runtime_error("kmeans: Lloyd iteration did not converge in " +
                               std::to_string(kMaxIterations) + " iterations");
    }
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

// GCBD (Grid-based Connectivity using Union-Find) -- WITHDRAWN from the comparison (author,
// 27.09.26): a plain 6-connected union-find CCL, identical in output to CC3D, and not the
// density-based grid clustering of Du & Wu (Entropy 24:1606, 2022) it was attributed to.
// Kept for the record and the exact-partition tests; not run in the campaigns.
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

namespace {

// ── HDBSCAN* ──────────────────────────────────────────────────────────────────
//
// As PUBLISHED: R.J.G.B. Campello, D. Moulavi, J. Sander, "Density-based clustering based on
// hierarchical density estimates", PAKDD 2013, LNCS 7819:160-172,
// doi:10.1007/978-3-642-37456-2_14. Audit 27.09.26 (decision D4/D5 revised): replaces a
// global cut at 1.5 x the median of an explicit N^2 edge list (finding S1-9), which was not
// HDBSCAN* and needed 16*N^2/2 bytes.
//
//   Definition 5  core distance d_core(x): distance to the m_pts-nearest neighbour, x
//                 included (m_pts = minSamples).
//   Definition 7  d_mreach(p, q) = max{d_core(p), d_core(q), d(p, q)}.
//   Algorithm 1   (2) an MST of the complete mutual-reachability graph, by "Prim's algorithm
//                 based on an ordinary list search (instead of a heap)" -- O(n^2) time, O(n)
//                 memory, no edge list; (3) extended with a self edge per object weighted by
//                 its core distance; (4) edges removed in decreasing order of weight, "in case
//                 of ties, edges must be removed simultaneously".
//   Algorithm 2   with m_clSize (= minClusterSize): after each removal every affected cluster
//                 is relabelled by its connected subcomponents; a subcomponent is spurious if
//                 it has fewer than m_clSize objects (or, for m_clSize = 1, no edge); all
//                 spurious -> the cluster disappears; one not spurious -> it shrinks, keeping
//                 its label; two or more -> a true split into new clusters.
//   Eq. (3)       S(C) = sum over x in C of (1/eps_min(x, C) - 1/eps_max(C)), eps_max(C) the
//                 level at which C appears (infinity for the root), eps_min(x, C) the level
//                 at which x leaves C (becomes noise, or C splits or disappears).
//   Algorithm 3   optimal flat selection over all clusters EXCEPT the root: bottom-up,
//                 S_hat(C) = S(C) for a leaf, else max{S(C), sum of the children's S_hat},
//                 keeping C when S(C) >= that sum; then the shallowest kept cluster per branch.
// A selected cluster's objects are those it held when it appeared; all other objects are
// noise and are not counted. Distances are Euclidean between voxel centres in voxel units,
// handled as exact integer squares so that ties are exact (they are frequent on a grid).
// m_pts must be >= 2 (with m_pts = 1 every core distance is 0 and every stability infinite).
// Under PBC xyz distances are minimum-image.
template <bool Periodic>
ClusterResult hdbscanImpl(
    int nx, int ny, int nz,
    int minClusterSize, int minSamples,
    const std::vector<uint8_t>& occupancy,
    std::vector<uint8_t>& visited) {

  ScopedTimer timer;
  ClusterResult result;

  if (minSamples < 2 || minClusterSize < 1) {
    throw std::invalid_argument("hdbscan: needs minSamples (m_pts) >= 2 and minClusterSize >= 1");
  }
  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  std::vector<int> occVox;
  for (std::size_t i = 0; i < totalSize; ++i)
    if (occupancy[i] == 1) occVox.push_back(static_cast<int>(i));
  const int n = static_cast<int>(occVox.size());
  if (n < minSamples) {  // no object has an m_pts-nearest neighbour: everything is noise
    result.elapsedMs = timer.elapsedMilliseconds();
    return result;
  }
  const int planeYZ = ny * nz;
  std::vector<int> px(static_cast<std::size_t>(n)), py(static_cast<std::size_t>(n)), pz(static_cast<std::size_t>(n));
  for (int c = 0; c < n; ++c) {
    const int v = occVox[static_cast<std::size_t>(c)];
    px[static_cast<std::size_t>(c)] = v / planeYZ;
    py[static_cast<std::size_t>(c)] = (v % planeYZ) / nz;
    pz[static_cast<std::size_t>(c)] = v % nz;
  }
  auto axis = [](int d, int len) {
    if constexpr (Periodic) {
      d = ((d % len) + len) % len;
      return std::min(d, len - d);
    } else {
      return d < 0 ? -d : d;
    }
  };
  auto dist2 = [&](int a, int b) -> long long {
    const long long dx = axis(px[static_cast<std::size_t>(a)] - px[static_cast<std::size_t>(b)], nx);
    const long long dy = axis(py[static_cast<std::size_t>(a)] - py[static_cast<std::size_t>(b)], ny);
    const long long dz = axis(pz[static_cast<std::size_t>(a)] - pz[static_cast<std::size_t>(b)], nz);
    return dx * dx + dy * dy + dz * dz;
  };

  // --- Definition 5: core distances (squared), by Chebyshev shells of growing radius s. After
  // shell s every occupied voxel within Chebyshev (min-image) distance s has been seen, and any
  // other lies farther than s, so the m-th smallest seen value is exact once it is <= s^2.
  std::vector<long long> core2(static_cast<std::size_t>(n));
  {
    const int maxShell = Periodic ? std::max({nx, ny, nz}) / 2 + 1 : std::max({nx, ny, nz});
    std::vector<long long> seen;
    std::vector<int> seenVox;
    for (int c = 0; c < n; ++c) {
      seen.clear();
      seenVox.clear();
      const int x = px[static_cast<std::size_t>(c)], y = py[static_cast<std::size_t>(c)], z = pz[static_cast<std::size_t>(c)];
      long long kth = std::numeric_limits<long long>::max();
      for (int s = 0; s <= maxShell; ++s) {
        for (int di = -s; di <= s; ++di)
          for (int dj = -s; dj <= s; ++dj)
            for (int dk = -s; dk <= s; ++dk) {
              if (std::max({std::abs(di), std::abs(dj), std::abs(dk)}) != s) continue;
              int a = x + di, b = y + dj, d = z + dk;
              if constexpr (Periodic) {
                a = ((a % nx) + nx) % nx;
                b = ((b % ny) + ny) % ny;
                d = ((d % nz) + nz) % nz;
              } else {
                if (a < 0 || a >= nx || b < 0 || b >= ny || d < 0 || d >= nz) continue;
              }
              const int w = static_cast<int>(idx3(a, b, d, ny, nz));
              if (occupancy[static_cast<std::size_t>(w)] != 1) continue;
              if constexpr (Periodic) {  // a wrapped offset can revisit a voxel in a small cell
                if (std::find(seenVox.begin(), seenVox.end(), w) != seenVox.end()) continue;
                seenVox.push_back(w);
              }
              const long long dx = axis(a - x, nx), dy = axis(b - y, ny), dz = axis(d - z, nz);
              seen.push_back(dx * dx + dy * dy + dz * dz);
            }
        if (static_cast<int>(seen.size()) >= minSamples) {
          std::nth_element(seen.begin(), seen.begin() + (minSamples - 1), seen.end());
          kth = seen[static_cast<std::size_t>(minSamples - 1)];
          if (kth <= static_cast<long long>(s) * s) break;
        }
      }
      core2[static_cast<std::size_t>(c)] = kth;
    }
  }

  // --- Algorithm 1 step 2: MST of the mutual-reachability graph, Prim with a list search.
  struct TreeEdge { int u, v; long long w; };
  std::vector<TreeEdge> mst;
  mst.reserve(static_cast<std::size_t>(n) - 1);
  {
    std::vector<long long> key(static_cast<std::size_t>(n), std::numeric_limits<long long>::max());
    std::vector<int> from(static_cast<std::size_t>(n), -1);
    std::vector<char> inTree(static_cast<std::size_t>(n), 0);
    int u = 0;
    for (int added = 0; added < n; ++added) {
      inTree[static_cast<std::size_t>(u)] = 1;
      if (from[static_cast<std::size_t>(u)] >= 0) mst.push_back({from[static_cast<std::size_t>(u)], u, key[static_cast<std::size_t>(u)]});
      int next = -1;
      long long best = std::numeric_limits<long long>::max();
      const long long cu = core2[static_cast<std::size_t>(u)];
      for (int v = 0; v < n; ++v) {
        if (inTree[static_cast<std::size_t>(v)]) continue;
        const long long w = std::max({cu, core2[static_cast<std::size_t>(v)], dist2(u, v)});
        if (w < key[static_cast<std::size_t>(v)]) { key[static_cast<std::size_t>(v)] = w; from[static_cast<std::size_t>(v)] = u; }
        if (key[static_cast<std::size_t>(v)] < best) { best = key[static_cast<std::size_t>(v)]; next = v; }
      }
      if (next < 0) break;
      u = next;
    }
  }

  // --- Algorithm 1 steps 3-4 with Algorithm 2: top-down removal of MST_ext edges by level.
  std::vector<std::vector<std::pair<int, long long>>> adj(static_cast<std::size_t>(n));
  for (const auto& e : mst) {
    adj[static_cast<std::size_t>(e.u)].push_back({e.v, e.w});
    adj[static_cast<std::size_t>(e.v)].push_back({e.u, e.w});
  }
  std::vector<long long> levels;
  levels.reserve(mst.size() + static_cast<std::size_t>(n));
  for (const auto& e : mst) levels.push_back(e.w);
  for (long long c2 : core2) levels.push_back(c2);
  std::sort(levels.begin(), levels.end(), std::greater<long long>());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  std::vector<int> edgeOrder(mst.size()), vertexOrder(static_cast<std::size_t>(n));
  std::iota(edgeOrder.begin(), edgeOrder.end(), 0);
  std::iota(vertexOrder.begin(), vertexOrder.end(), 0);
  std::sort(edgeOrder.begin(), edgeOrder.end(), [&](int a, int b) { return mst[static_cast<std::size_t>(a)].w > mst[static_cast<std::size_t>(b)].w; });
  std::sort(vertexOrder.begin(), vertexOrder.end(), [&](int a, int b) { return core2[static_cast<std::size_t>(a)] > core2[static_cast<std::size_t>(b)]; });

  struct Cluster {
    int parent;
    double lambdaBirth;  // 1 / eps_max(C); 0 for the root
    double stability{0.0};
    std::vector<int> children;
    std::vector<int> birthMembers;
  };
  std::vector<Cluster> cl;
  std::vector<std::vector<int>> members;  // current members of each cluster
  std::vector<int> clusterOf(static_cast<std::size_t>(n), 0);
  cl.push_back({-1, 0.0, 0.0, {}, {}});
  members.emplace_back(static_cast<std::size_t>(n));
  std::iota(members[0].begin(), members[0].end(), 0);

  std::vector<int> comp(static_cast<std::size_t>(n), -1), stack, affected;
  std::size_t ei = 0, vi = 0;
  for (const long long w : levels) {
    const double lambda = 1.0 / std::sqrt(static_cast<double>(w));
    affected.clear();
    for (; ei < edgeOrder.size() && mst[static_cast<std::size_t>(edgeOrder[ei])].w == w; ++ei) {
      const int c = clusterOf[static_cast<std::size_t>(mst[static_cast<std::size_t>(edgeOrder[ei])].u)];
      if (c >= 0) affected.push_back(c);
    }
    for (; vi < vertexOrder.size() && core2[static_cast<std::size_t>(vertexOrder[vi])] == w; ++vi) {
      const int c = clusterOf[static_cast<std::size_t>(vertexOrder[vi])];
      if (c >= 0) affected.push_back(c);
    }
    std::sort(affected.begin(), affected.end());
    affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
    for (const int c : affected) {
      // Subcomponents of cluster c under the remaining edges (weight < w).
      std::vector<std::vector<int>> subs;
      for (const int v : members[static_cast<std::size_t>(c)]) comp[static_cast<std::size_t>(v)] = -1;
      for (const int v0 : members[static_cast<std::size_t>(c)]) {
        if (comp[static_cast<std::size_t>(v0)] >= 0) continue;
        const int id = static_cast<int>(subs.size());
        subs.emplace_back();
        stack.assign(1, v0);
        comp[static_cast<std::size_t>(v0)] = id;
        while (!stack.empty()) {
          const int v = stack.back();
          stack.pop_back();
          subs.back().push_back(v);
          for (const auto& [q, ew] : adj[static_cast<std::size_t>(v)]) {
            if (ew >= w || clusterOf[static_cast<std::size_t>(q)] != c || comp[static_cast<std::size_t>(q)] >= 0) continue;
            comp[static_cast<std::size_t>(q)] = id;
            stack.push_back(q);
          }
        }
      }
      std::vector<int> keep;
      for (int s = 0; s < static_cast<int>(subs.size()); ++s) {
        const auto& sub = subs[static_cast<std::size_t>(s)];
        const bool hasEdge = sub.size() >= 2 || core2[static_cast<std::size_t>(sub[0])] < w;
        const bool spurious = static_cast<int>(sub.size()) < minClusterSize || (minClusterSize == 1 && !hasEdge);
        if (!spurious) keep.push_back(s);
      }
      Cluster& C = cl[static_cast<std::size_t>(c)];
      const double gain = lambda - C.lambdaBirth;  // 1/eps_min - 1/eps_max for an object leaving now
      if (keep.size() == 1) {                        // shrink: spurious parts leave as noise
        std::vector<int> stay = subs[static_cast<std::size_t>(keep[0])];
        const std::size_t left = members[static_cast<std::size_t>(c)].size() - stay.size();
        cl[static_cast<std::size_t>(c)].stability += static_cast<double>(left) * gain;
        for (const int v : members[static_cast<std::size_t>(c)]) clusterOf[static_cast<std::size_t>(v)] = -1;
        for (const int v : stay) clusterOf[static_cast<std::size_t>(v)] = c;
        members[static_cast<std::size_t>(c)].swap(stay);
      } else {                                        // disappears (0) or true split (>= 2)
        cl[static_cast<std::size_t>(c)].stability += static_cast<double>(members[static_cast<std::size_t>(c)].size()) * gain;
        for (const int v : members[static_cast<std::size_t>(c)]) clusterOf[static_cast<std::size_t>(v)] = -1;
        members[static_cast<std::size_t>(c)].clear();
        for (const int s : keep) {
          const int id = static_cast<int>(cl.size());
          cl.push_back({c, lambda, 0.0, {}, subs[static_cast<std::size_t>(s)]});
          cl[static_cast<std::size_t>(c)].children.push_back(id);
          for (const int v : subs[static_cast<std::size_t>(s)]) clusterOf[static_cast<std::size_t>(v)] = id;
          members.push_back(subs[static_cast<std::size_t>(s)]);
        }
      }
    }
  }

  // --- Algorithm 3: optimal flat selection, root excluded.
  const int nc = static_cast<int>(cl.size());
  std::vector<double> sHat(static_cast<std::size_t>(nc), 0.0);
  std::vector<char> keepC(static_cast<std::size_t>(nc), 0);
  for (int c = nc - 1; c >= 1; --c) {  // children are created after their parents
    const Cluster& C = cl[static_cast<std::size_t>(c)];
    if (C.children.empty()) { sHat[static_cast<std::size_t>(c)] = C.stability; keepC[static_cast<std::size_t>(c)] = 1; continue; }
    double sum = 0.0;
    for (const int ch : C.children) sum += sHat[static_cast<std::size_t>(ch)];
    if (C.stability < sum) { sHat[static_cast<std::size_t>(c)] = sum; keepC[static_cast<std::size_t>(c)] = 0; }
    else { sHat[static_cast<std::size_t>(c)] = C.stability; keepC[static_cast<std::size_t>(c)] = 1; }
  }
  std::vector<int> selected, todo(cl[0].children.rbegin(), cl[0].children.rend());
  while (!todo.empty()) {
    const int c = todo.back();
    todo.pop_back();
    if (keepC[static_cast<std::size_t>(c)]) { selected.push_back(c); continue; }
    for (auto it = cl[static_cast<std::size_t>(c)].children.rbegin(); it != cl[static_cast<std::size_t>(c)].children.rend(); ++it) todo.push_back(*it);
  }

  for (const int c : selected) {
    const auto& mem = cl[static_cast<std::size_t>(c)].birthMembers;
    const int sz = static_cast<int>(mem.size());
    for (const int v : mem) visited[static_cast<std::size_t>(occVox[static_cast<std::size_t>(v)])] = 1;
    result.visitedVoxels += static_cast<std::size_t>(sz);
    result.nclusters++;
    result.clusterSizes.push_back(sz);
    result.maxCluster = std::max(result.maxCluster, sz);
  }
  std::sort(result.clusterSizes.begin(), result.clusterSizes.end(), std::greater<int>());
  result.elapsedMs = timer.elapsedMilliseconds();
  return result;
}

}  // namespace

ClusterResult hdbscan(int nx, int ny, int nz, int minClusterSize, int minSamples,
                      const std::vector<uint8_t>& occupancy, std::vector<uint8_t>& visited,
                      bool periodic) {
  return periodic ? hdbscanImpl<true>(nx, ny, nz, minClusterSize, minSamples, occupancy, visited)
                  : hdbscanImpl<false>(nx, ny, nz, minClusterSize, minSamples, occupancy, visited);
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

// ── CC3D, optimized track: SAUF ──────────────────────────────────────────────
//
// Scan plus Array-based Union-Find (SAUF) as PUBLISHED: K. Wu, E. Otoo, A. Shoshani,
// "Optimizing connected component labeling algorithms", SPIE 5747 (2005), App. A, and
// K. Wu, E. Otoo, K. Suzuki, "Optimizing two-pass connected-component labeling
// algorithms", Pattern Anal. Appl. 12(2):117-135, doi:10.1007/s10044-008-0109-y, §3.3.
// Audit 27.09.26: this replaces a variant with a rolling two-plane label buffer and an
// occupied-voxel-only labeling pass, neither of which is in the sources.
//
//   scanning phase  raster order over every voxel; the scan mask is the voxel's already
//                   scanned neighbours (3 for 6-connectivity, 13 for 26). Eq. (8): no
//                   labelled neighbour -> new label l (P[l] <- l); otherwise
//                   L[e] <- min_i findRoot(P, L[i]), then setRoot(P, L[i], L[e]) for all i.
//                   This is the scan S0; the decision tree (S1) of the papers is for 2-D
//                   8-connectivity and has no 3-D form in them.
//   union-find      the papers' array: findRoot follows P[i] < i to the root; setRoot
//                   points every node of a path at the root; union keeps the smaller root
//                   and applies setRoot to both nodes (so P[i] <= i always).
//   analysis phase  flattenL: one pass over P assigning consecutive final labels.
//   labeling phase  a second pass over the whole label image, L[e] <- P[L[e]].
// The provisional label image L is grid-sized, as in the papers. Label 0 is background.
// Under PBC xyz, the neighbour pairs that wrap across a face are united (the papers'
// union) after the scan and before the analysis phase -- an addition for periodic grids.
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

  std::vector<unsigned> L(totalSize, 0u);  // provisional label image
  std::vector<unsigned> P;                 // equivalence array; P[0] = background
  P.reserve(1024);
  P.push_back(0u);

  auto findRoot = [&P](unsigned i) {
    while (P[i] < i) i = P[i];
    return i;
  };
  auto setRoot = [&P](unsigned i, unsigned root) {
    while (P[i] < i) {
      const unsigned j = P[i];
      P[i] = root;
      i = j;
    }
    P[i] = root;
  };
  auto unite = [&](unsigned i, unsigned j) {
    unsigned root = findRoot(i);
    if (i != j) {
      const unsigned rootj = findRoot(j);
      if (root > rootj) root = rootj;
      setRoot(j, root);
    }
    setRoot(i, root);
    return root;
  };

  // Scan mask: the neighbours already visited in (i,j,k) raster order.
  constexpr int back6[3][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
  constexpr int back26[13][3] = {
      {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
      {-1, -1, 0}, {-1, 1, 0}, {-1, 0, -1}, {-1, 0, 1},
      {0, -1, -1}, {0, -1, 1},
      {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}, {-1, 1, 1}};
  const int (*back)[3] = (connectivity == 26) ? back26 : back6;
  const int numBack = (connectivity == 26) ? 13 : 3;

  // --- scanning phase (eq. 8) --------------------------------------------------------
  unsigned nbrLabel[13];
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t e = idx3(i, j, k, ny, nz);
        if (occupancy[e] != 1) continue;
        int n = 0;
        for (int d = 0; d < numBack; ++d) {
          const int ni = i + back[d][0], nj = j + back[d][1], nk = k + back[d][2];
          if (ni < 0 || nj < 0 || nj >= ny || nk < 0 || nk >= nz) continue;
          const unsigned l = L[idx3(ni, nj, nk, ny, nz)];
          if (l != 0u) nbrLabel[n++] = l;
        }
        if (n == 0) {
          const unsigned l = static_cast<unsigned>(P.size());
          P.push_back(l);
          L[e] = l;
          continue;
        }
        unsigned m = findRoot(nbrLabel[0]);
        for (int t = 1; t < n; ++t) m = std::min(m, findRoot(nbrLabel[t]));
        for (int t = 0; t < n; ++t) setRoot(nbrLabel[t], m);
        L[e] = m;
      }
    }
  }

  if constexpr (Periodic) {
    // Pairs across a face: from every voxel on a face, every stencil neighbour that lies
    // outside the grid, wrapped. In-grid pairs were all seen by the scan.
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < nz; ++k) {
          if (i != 0 && i != nx - 1 && j != 0 && j != ny - 1 && k != 0 && k != nz - 1) continue;
          const unsigned a = L[idx3(i, j, k, ny, nz)];
          if (a == 0u) continue;
          for (int di = -1; di <= 1; ++di)
            for (int dj = -1; dj <= 1; ++dj)
              for (int dk = -1; dk <= 1; ++dk) {
                const int m = std::abs(di) + std::abs(dj) + std::abs(dk);
                if (m == 0 || (connectivity != 26 && m != 1)) continue;
                int ni = i + di, nj = j + dj, nk = k + dk;
                if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) continue;
                ni = (ni + nx) % nx;
                nj = (nj + ny) % ny;
                nk = (nk + nz) % nz;
                const unsigned b = L[idx3(ni, nj, nk, ny, nz)];
                if (b != 0u) unite(a, b);
              }
        }
      }
    }
  }

  // --- analysis phase: flattenL (consecutive final labels 1..count) ------------------
  unsigned count = 0;
  for (std::size_t l = 1; l < P.size(); ++l) {
    if (P[l] < l) P[l] = P[P[l]];
    else P[l] = ++count;
  }

  // --- labeling phase: second pass over the image -------------------------------------
  std::vector<int> sizes(count, 0);
  for (std::size_t e = 0; e < totalSize; ++e) {
    if (L[e] == 0u) continue;
    const unsigned f = P[L[e]];
    L[e] = f;
    ++sizes[f - 1];
    if (labels) (*labels)[e] = static_cast<int>(f - 1);
    visited[e] = 1;
    result.visitedVoxels++;
  }

  result.nclusters = static_cast<int>(count);
  for (int s : sizes) {
    result.clusterSizes.push_back(s);
    result.maxCluster = std::max(result.maxCluster, s);
  }
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



// ── RLE-CCL, optimized track: He, Chao & Suzuki 2008 ─────────────────────────
//
// The run-based two-scan labeling algorithm as PUBLISHED: L. He, Y. Chao, K. Suzuki,
// "A run-based two-scan labeling algorithm", IEEE Trans. Image Process. 17(5):749-756,
// 2008, doi:10.1109/TIP.2008.919369, §II-III, generalised from 2-D rows to 3-D rows.
// Audit 27.09.26: replaces a union-find over run ids with a second pass over the runs,
// which is not the paper's method.
//
//   first scan    raster order. A run is a maximal block of occupied voxels along k in
//                 row (i,j). Runs connected to the current run are those of the already
//                 scanned neighbouring rows -- (i-1,j) and (i,j-1) in 3-D -- whose k-range
//                 overlaps it (6-connectivity; the paper's 2-D case is 8-connectivity).
//                 No connected run: a new provisional label l and a new provisional label
//                 list, rtable[l] = l, next[l] = -1, tail[l] = l. Otherwise the run takes
//                 the provisional label of the first connected run in scan order (the
//                 paper's "upper left-most run") and resolve() merges the lists of every
//                 other connected run into it. Every voxel of the run is written with its
//                 provisional label (label image b).
//   label lists   the paper's connection lists: rtable (representative = smallest label =
//                 list head), next, tail. merge(u, v), u < v: every label of list v gets
//                 rtable = u, next[tail[u]] = v, tail[u] = tail[v]. resolve(x, y): u =
//                 rtable[x], v = rtable[y]; merge the larger into the smaller; nothing if
//                 equal.
//   second scan   over the whole image, b <- rtable[b].
// The runs of the previous rows are kept per row (the paper keeps them in a start/end queue,
// which serves the single previous row of 2-D). Label 0 is background. Under PBC xyz the
// voxel pairs that wrap across a face are resolved after the first scan (addition).
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

  std::vector<unsigned> b(totalSize, 0u);  // label image
  // Provisional label lists (index 0 = background, rtable[0] = 0).
  std::vector<unsigned> rtable{0u}, tail{0u};
  std::vector<long long> next{-1};

  auto merge = [&](unsigned u, unsigned v) {  // u < v: list v into list u
    for (long long i = v; i != -1; i = next[static_cast<std::size_t>(i)]) rtable[static_cast<std::size_t>(i)] = u;
    next[tail[u]] = v;
    tail[u] = tail[v];
  };
  auto resolve = [&](unsigned x, unsigned y) {
    const unsigned u = rtable[x], v = rtable[y];
    if (u < v) merge(u, v);
    else if (u > v) merge(v, u);
  };

  struct Run { int kStart, kEnd; unsigned label; };
  std::vector<std::vector<Run>> prevSlice(static_cast<std::size_t>(ny)), currSlice(static_cast<std::size_t>(ny));
  std::vector<unsigned> connected;

  // --- first scan ----------------------------------------------------------------------
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      auto& row = currSlice[static_cast<std::size_t>(j)];
      row.clear();
      const std::size_t rowBase = idx3(i, j, 0, ny, nz);
      int k = 0;
      while (k < nz) {
        if (occupancy[rowBase + static_cast<std::size_t>(k)] != 1) { ++k; continue; }
        const int ks = k;
        while (k < nz && occupancy[rowBase + static_cast<std::size_t>(k)] == 1) ++k;
        const int ke = k - 1;
        // Connected runs of the scanned neighbour rows, in scan order: (i-1,j) then (i,j-1).
        connected.clear();
        if (i > 0)
          for (const Run& r : prevSlice[static_cast<std::size_t>(j)])
            if (r.kEnd >= ks && r.kStart <= ke) connected.push_back(r.label);
        if (j > 0)
          for (const Run& r : currSlice[static_cast<std::size_t>(j - 1)])
            if (r.kEnd >= ks && r.kStart <= ke) connected.push_back(r.label);
        unsigned lab;
        if (connected.empty()) {
          lab = static_cast<unsigned>(rtable.size());
          rtable.push_back(lab);
          next.push_back(-1);
          tail.push_back(lab);
        } else {
          lab = connected[0];
          for (std::size_t t = 1; t < connected.size(); ++t) resolve(lab, connected[t]);
        }
        for (int kk = ks; kk <= ke; ++kk) b[rowBase + static_cast<std::size_t>(kk)] = lab;
        row.push_back(Run{ks, ke, lab});
      }
    }
    std::swap(prevSlice, currSlice);
  }

  if constexpr (Periodic) {
    // Voxel pairs across a face (6-connectivity): x, y and z faces.
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
        const unsigned a = b[idx3(0, j, k, ny, nz)], c = b[idx3(nx - 1, j, k, ny, nz)];
        if (a && c) resolve(a, c);
      }
    for (int i = 0; i < nx; ++i)
      for (int k = 0; k < nz; ++k) {
        const unsigned a = b[idx3(i, 0, k, ny, nz)], c = b[idx3(i, ny - 1, k, ny, nz)];
        if (a && c) resolve(a, c);
      }
    for (int i = 0; i < nx; ++i)
      for (int j = 0; j < ny; ++j) {
        const unsigned a = b[idx3(i, j, 0, ny, nz)], c = b[idx3(i, j, nz - 1, ny, nz)];
        if (a && c) resolve(a, c);
      }
  }

  // --- second scan: b <- rtable[b] over the whole image --------------------------------
  // Output contract (dense ids, sizes) is tallied in the same pass: representatives are
  // numbered in order of first appearance, as every method's labels are.
  std::vector<int> dense(rtable.size(), -1);
  std::vector<int> sizes;
  for (std::size_t e = 0; e < totalSize; ++e) {
    if (b[e] == 0u) continue;
    const unsigned rep = rtable[b[e]];
    b[e] = rep;
    int& id = dense[rep];
    if (id < 0) { id = static_cast<int>(sizes.size()); sizes.push_back(0); }
    ++sizes[static_cast<std::size_t>(id)];
    if (labels) (*labels)[e] = id;
    visited[e] = 1;
    result.visitedVoxels++;
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
// Voxel Cloud Connectivity Segmentation as PUBLISHED: J. Papon, A. Abramov, M. Schoeler,
// F. Wörgötter, "Voxel Cloud Connectivity Segmentation - Supervoxels for Point Clouds",
// CVPR 2013, doi:10.1109/CVPR.2013.264, §3.1-3.4 (audit decision D6 revised: the method
// follows the paper, not PCL's SupervoxelClustering code, which departs from it). The
// voxel grid is the paper's voxel cloud with R_voxel = 1 voxel; R_seed = S voxels.
//
//   §3.1 adjacency   26-adjacency ("whenever we refer to adjacent voxels, we are speaking
//                    of 26-adjacency").
//   §3.2 seeding     one candidate per occupied seed voxel (S^3 cell): "the voxel in the
//                    cloud nearest to the center of each occupied seeding voxel".
//        filtering   delete candidates that do "not have at least as many voxels as would be
//                    occupied by a planar surface intersecting with half of the search volume":
//                    kept iff #{occupied v : |v - c| <= R_search} >= pi R_search^2 / R_voxel^2.
//                    The paper gives no value for R_search; 0.5 R_seed is taken from PCL (by
//                    the paper's first author), as is area/voxel-area for "voxels occupied by
//                    a planar slice" -- without PCL's extra factor 0.05, which PCL itself
//                    notes is "smaller than the value mentioned in the original paper".
//        seed shift  "to the connected voxel within the search volume which has the smallest
//                    gradient", the gradient being a CIELab colour gradient (eq. 2). On a
//                    binary occupancy grid there is no colour, the gradient is 0 everywhere
//                    and the seed stays where it is.
//   §3.3 centre      initialised as "the center ... of the seed voxel and connected neighbors
//                    within 2 voxels" (read: occupied voxels within 2 steps of the seed in the
//                    adjacency graph). Distance D (eq. 4) keeps only its spatial term: colour
//                    (D_c) and FPFH (D_HiK) need colour and surface normals, which a binary
//                    occupancy grid does not have; D is then proportional to |v - centre|.
//        search      limited "so that it ends at the neighboring cluster centers", the farthest
//        volume      point lying at sqrt(3) R_seed: the cube |v - centre|_k <= R_seed per axis.
//   §3.4 flow        per iteration, per supervoxel, a breadth-first search of the adjacency
//                    graph from "the voxel nearest the cluster center", one level at a time for
//                    ALL supervoxels before going deeper. A voxel reached takes the label "if the
//                    distance is the smallest this voxel has seen"; only then are its neighbours
//                    "which are further from the center" added to this supervoxel's queue,
//                    never twice in one iteration. A supervoxel's search ends when it runs out
//                    of voxels or "none of the nodes searched in the current level were set to
//                    its label". Then every centre becomes "the mean of all its constituents".
//        iterations  "until the cluster centers stabilize, or for a fixed number of iterations
//                    ... five iterations for all presented results": at most 5, stopping early
//                    when no centre moves (the next iteration would repeat this one).
//
// Readings stated because the paper does not spell them out: "smallest this voxel has seen"
// restarts every iteration (the local k-means of SLIC [1], which §3.4 builds on); ties keep
// the earlier supervoxel; supervoxels and queues are processed in voxel-index order.
// Voxels no search reaches stay unlabelled in the paper; here each 26-connected remainder
// becomes one cluster so that every method is tallied over the same occupied set.
// Under PBC xyz adjacency wraps and every distance is a minimum image; the seed-cell grid is
// not periodic (a cell cut by a face is a smaller cell).
//
// Implementation: no grid-sized array; the occupied voxels are listed once in raster order
// and all per-voxel state is indexed by position in that list (voxel -> position by binary
// search).
template <bool Periodic>
ClusterResult vccsOptimizedImpl(int nx, int ny, int nz, double seedResolution,
                                const std::vector<uint8_t>& occupancy,
                                std::vector<uint8_t>& visited) {
  ScopedTimer timer;
  ClusterResult result;

  const std::size_t totalSize = static_cast<std::size_t>(nx) * ny * nz;
  std::fill(visited.begin(), visited.end(), 0);

  const int S = std::max(1, static_cast<int>(std::round(seedResolution)));

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
  std::vector<std::array<int, 3>> pos(static_cast<std::size_t>(nOcc));
  for (int c = 0; c < nOcc; ++c) {
    const int v = occVox[static_cast<std::size_t>(c)];
    const int i = v / planeYZ, rem = v - i * planeYZ;
    pos[static_cast<std::size_t>(c)] = {i, rem / nz, rem % nz};
  }
  // Occupied voxel at (i,j,k) -> position in the list, -1 if empty or (non-periodic) outside.
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
  const int dims[3] = {nx, ny, nz};
  auto delta = [](double a, double b, int n) {
    double d = a - b;
    if constexpr (Periodic) d -= n * std::round(d / n);
    return d;
  };
  using P3 = std::array<double, 3>;
  auto dist2 = [&](int c, const P3& p) {
    const auto& q = pos[static_cast<std::size_t>(c)];
    double s = 0.0;
    for (int a = 0; a < 3; ++a) {
      const double d = delta(q[a], p[a], dims[a]);
      s += d * d;
    }
    return s;
  };
  std::vector<std::array<int, 3>> adj26;
  for (int di = -1; di <= 1; ++di)
    for (int dj = -1; dj <= 1; ++dj)
      for (int dk = -1; dk <= 1; ++dk)
        if (di || dj || dk) adj26.push_back({di, dj, dk});
  // Nearest occupied voxel to p: an initial guess bounds the cube that must be searched.
  auto nearestTo = [&](const P3& p, int guess) {
    int best = guess;
    double bestD2 = dist2(guess, p);
    const int reach = static_cast<int>(std::ceil(std::sqrt(bestD2))) + 1;
    const int ci = static_cast<int>(std::floor(p[0])), cj = static_cast<int>(std::floor(p[1])),
              ck = static_cast<int>(std::floor(p[2]));
    for (int i = ci - reach; i <= ci + reach + 1; ++i)
      for (int j = cj - reach; j <= cj + reach + 1; ++j)
        for (int k = ck - reach; k <= ck + reach + 1; ++k) {
          const int c = occupiedAt(i, j, k);
          if (c < 0 || c == best) continue;
          const double d2 = dist2(c, p);
          if (d2 < bestD2 || (d2 == bestD2 && c < best)) { best = c; bestD2 = d2; }
        }
    return best;
  };

  // --- §3.2 seeding: nearest occupied voxel to each occupied seed-cell centre -------------
  std::vector<std::pair<std::array<int, 3>, int>> byCell;
  byCell.reserve(static_cast<std::size_t>(nOcc));
  for (int c = 0; c < nOcc; ++c) {
    const auto& q = pos[static_cast<std::size_t>(c)];
    byCell.push_back({{q[0] / S, q[1] / S, q[2] / S}, c});
  }
  std::sort(byCell.begin(), byCell.end());
  std::vector<int> candidates;
  const double half = 0.5 * (S - 1);
  for (std::size_t g = 0; g < byCell.size();) {
    const auto cell = byCell[g].first;
    const P3 centre{cell[0] * S + half, cell[1] * S + half, cell[2] * S + half};
    const int first = byCell[g].second;
    while (g < byCell.size() && byCell[g].first == cell) ++g;
    candidates.push_back(nearestTo(centre, first));
  }
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

  // --- §3.2 filtering ---------------------------------------------------------------------
  const double rSearch = 0.5 * S;
  const double minVoxels = 3.14159265358979323846 * rSearch * rSearch;
  const int reachR = static_cast<int>(std::floor(rSearch));
  std::vector<int> seeds, within;
  for (int c : candidates) {
    const auto& q = pos[static_cast<std::size_t>(c)];
    const P3 qc{double(q[0]), double(q[1]), double(q[2])};
    within.clear();
    for (int di = -reachR; di <= reachR; ++di)
      for (int dj = -reachR; dj <= reachR; ++dj)
        for (int dk = -reachR; dk <= reachR; ++dk) {
          const int u = occupiedAt(q[0] + di, q[1] + dj, q[2] + dk);
          if (u >= 0 && dist2(u, qc) <= rSearch * rSearch) within.push_back(u);
        }
    // Distinct voxels: in a periodic axis shorter than 2 R_search two offsets reach one voxel.
    std::sort(within.begin(), within.end());
    const auto count = std::unique(within.begin(), within.end()) - within.begin();
    if (count >= minVoxels) seeds.push_back(c);
  }
  result.seedCandidates = static_cast<int>(candidates.size());
  result.seedsPlaced = static_cast<int>(seeds.size());
  result.seedPruneThreshold = minVoxels;

  // --- §3.3 initial centres: seed and connected neighbours within 2 steps -----------------
  const int nSv = static_cast<int>(seeds.size());
  std::vector<P3> centre(static_cast<std::size_t>(nSv));
  std::vector<char> alive(static_cast<std::size_t>(nSv), 1);
  std::vector<int> nearGuess(seeds);  // a voxel near each centre, to start nearestTo()
  {
    std::vector<int> mark(static_cast<std::size_t>(nOcc), -1), layer, next, members;
    for (int s = 0; s < nSv; ++s) {
      const int c0 = seeds[static_cast<std::size_t>(s)];
      const auto& q0 = pos[static_cast<std::size_t>(c0)];
      members.assign(1, c0);
      layer = {c0};
      mark[static_cast<std::size_t>(c0)] = s;
      for (int step = 0; step < 2 && !layer.empty(); ++step) {
        next.clear();
        for (int c : layer) {
          const auto& q = pos[static_cast<std::size_t>(c)];
          for (const auto& o : adj26) {
            const int u = occupiedAt(q[0] + o[0], q[1] + o[1], q[2] + o[2]);
            if (u < 0 || mark[static_cast<std::size_t>(u)] == s) continue;
            mark[static_cast<std::size_t>(u)] = s;
            next.push_back(u);
            members.push_back(u);
          }
        }
        layer.swap(next);
      }
      // Summed in voxel order, like every centre update below, so a centre does not depend
      // on the order the neighbourhood happened to be discovered in.
      std::sort(members.begin(), members.end());
      P3 sum{0, 0, 0};
      for (int c : members) {
        const auto& q = pos[static_cast<std::size_t>(c)];
        for (int a = 0; a < 3; ++a) sum[a] += q0[a] + delta(q[a], q0[a], dims[a]);
      }
      const double n = static_cast<double>(members.size());
      centre[static_cast<std::size_t>(s)] = {sum[0] / n, sum[1] / n, sum[2] / n};
    }
  }

  // --- §3.4 flow-constrained clustering ---------------------------------------------------
  std::vector<int> label(static_cast<std::size_t>(nOcc), -1);
  std::vector<double> best(static_cast<std::size_t>(nOcc));
  std::vector<std::vector<int>> frontier(static_cast<std::size_t>(nSv));
  std::vector<std::unordered_set<int>> queued(static_cast<std::size_t>(nSv));
  std::vector<char> searching(static_cast<std::size_t>(nSv));
  std::vector<int> next;
  auto inSearchVolume = [&](int c, const P3& p) {
    const auto& q = pos[static_cast<std::size_t>(c)];
    for (int a = 0; a < 3; ++a)
      if (std::fabs(delta(q[a], p[a], dims[a])) > S) return false;
    return true;
  };
  for (int iter = 0; iter < 5; ++iter) {
    std::fill(label.begin(), label.end(), -1);
    std::fill(best.begin(), best.end(), std::numeric_limits<double>::max());
    for (int s = 0; s < nSv; ++s) {
      auto& F = frontier[static_cast<std::size_t>(s)];
      auto& Q = queued[static_cast<std::size_t>(s)];
      F.clear();
      Q.clear();
      searching[static_cast<std::size_t>(s)] = alive[static_cast<std::size_t>(s)];
      if (!alive[static_cast<std::size_t>(s)]) continue;
      const int start = nearestTo(centre[static_cast<std::size_t>(s)], nearGuess[static_cast<std::size_t>(s)]);
      F.push_back(start);
      Q.insert(start);
    }
    for (bool any = true; any;) {  // one level of every supervoxel's search per pass
      any = false;
      for (int s = 0; s < nSv; ++s) {
        if (!searching[static_cast<std::size_t>(s)]) continue;
        const P3& cen = centre[static_cast<std::size_t>(s)];
        auto& F = frontier[static_cast<std::size_t>(s)];
        auto& Q = queued[static_cast<std::size_t>(s)];
        next.clear();
        bool setAny = false;
        for (int v : F) {
          const double d = std::sqrt(dist2(v, cen));
          if (!(d < best[static_cast<std::size_t>(v)])) continue;
          best[static_cast<std::size_t>(v)] = d;
          label[static_cast<std::size_t>(v)] = s;
          setAny = true;
          const auto& q = pos[static_cast<std::size_t>(v)];
          for (const auto& o : adj26) {
            const int u = occupiedAt(q[0] + o[0], q[1] + o[1], q[2] + o[2]);
            if (u < 0 || Q.count(u) || !inSearchVolume(u, cen)) continue;
            if (!(std::sqrt(dist2(u, cen)) > d)) continue;  // only neighbours further out
            Q.insert(u);
            next.push_back(u);
          }
        }
        F.swap(next);
        if (!setAny || F.empty()) searching[static_cast<std::size_t>(s)] = 0;
        else any = true;
      }
    }
    // centres <- mean of constituents (minimum image about the old centre)
    std::vector<P3> sum(static_cast<std::size_t>(nSv), P3{0, 0, 0});
    std::vector<int> count(static_cast<std::size_t>(nSv), 0);
    for (int c = 0; c < nOcc; ++c) {
      const int s = label[static_cast<std::size_t>(c)];
      if (s < 0) continue;
      const P3& c0 = centre[static_cast<std::size_t>(s)];
      const auto& q = pos[static_cast<std::size_t>(c)];
      for (int a = 0; a < 3; ++a) sum[static_cast<std::size_t>(s)][a] += c0[a] + delta(q[a], c0[a], dims[a]);
      if (count[static_cast<std::size_t>(s)]++ == 0) nearGuess[static_cast<std::size_t>(s)] = c;
    }
    bool moved = false;
    for (int s = 0; s < nSv; ++s) {
      if (!alive[static_cast<std::size_t>(s)]) continue;
      const int n = count[static_cast<std::size_t>(s)];
      if (n == 0) { alive[static_cast<std::size_t>(s)] = 0; moved = true; continue; }
      const P3 c{sum[static_cast<std::size_t>(s)][0] / n, sum[static_cast<std::size_t>(s)][1] / n,
                 sum[static_cast<std::size_t>(s)][2] / n};
      if (c != centre[static_cast<std::size_t>(s)]) moved = true;
      centre[static_cast<std::size_t>(s)] = c;
    }
    if (!moved) break;
  }

  // --- tally; voxels no search reached: one cluster per 26-connected remainder -----------
  std::vector<int> assignment(static_cast<std::size_t>(nOcc), -1);
  std::vector<int> dense(static_cast<std::size_t>(nSv), -1);
  int clusterId = 0;
  for (int c = 0; c < nOcc; ++c) {
    const int s = label[static_cast<std::size_t>(c)];
    if (s < 0) continue;
    if (dense[static_cast<std::size_t>(s)] < 0) dense[static_cast<std::size_t>(s)] = clusterId++;
    assignment[static_cast<std::size_t>(c)] = dense[static_cast<std::size_t>(s)];
  }
  std::vector<int> stack;
  for (int c = 0; c < nOcc; ++c) {
    if (assignment[static_cast<std::size_t>(c)] != -1) continue;
    assignment[static_cast<std::size_t>(c)] = clusterId;
    stack.assign(1, c);
    while (!stack.empty()) {
      const int cur = stack.back();
      stack.pop_back();
      const auto& q = pos[static_cast<std::size_t>(cur)];
      for (const auto& o : adj26) {
        const int u = occupiedAt(q[0] + o[0], q[1] + o[1], q[2] + o[2]);
        if (u < 0 || assignment[static_cast<std::size_t>(u)] != -1) continue;
        assignment[static_cast<std::size_t>(u)] = clusterId;
        stack.push_back(u);
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
                            bool periodic) {
  return periodic ? vccsOptimizedImpl<true>(nx, ny, nz, seedResolution, occupancy, visited)
                  : vccsOptimizedImpl<false>(nx, ny, nz, seedResolution, occupancy, visited);
}

}  // namespace bls
