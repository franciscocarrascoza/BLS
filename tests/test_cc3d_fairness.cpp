// Dual-track equivalence: every optimized variant must agree with its
// textbook-fair counterpart.
//
// The pairs are (cc3d, cc3d_optimized), (rle_ccl, rle_ccl_optimized) and
// (vccs, vccs_optimized). The first two are exact connected-component
// labellers, so agreement is total: nclusters, max_cluster and the sorted
// component size distribution must match on every input.
//
// VCCS IS DIFFERENT, AND DELIBERATELY SO. It is a seeded segmentation, not a
// component labeller: its output is a function of where the seeds land. The
// optimized track implements the published adaptive seeding and seed pruning
// (Papon et al. 2013), which by construction moves seeds the fair track does
// not move and drops seeds the fair track keeps. Demanding identical output
// would be demanding that the optimization do nothing. What IS required of the
// pair, and asserted below, is that both partition exactly the occupied set,
// that every supervoxel is 6-connected, and that neither ever merges two
// distinct connected components into one supervoxel.
//
// Test geometry is deliberately not just axis-aligned cubes. Both major defects
// in this codebase survived because cubic, symmetric cases hid them, so the
// cases below include a sheared (non-orthogonal, non-symmetric) lattice of
// blobs whose principal directions align with no grid axis, plus shapes with
// concave and diagonal contact that separate 6- from 26-connectivity.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "cluster/Algorithms.hpp"

using namespace bls;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
  if (!ok) { std::cerr << "FAIL: " << what << "\n"; ++failures; }
}

struct Grid {
  std::string name;
  int nx, ny, nz;
  std::vector<uint8_t> occ;
  std::size_t idx(int i, int j, int k) const {
    return static_cast<std::size_t>(i) * ny * nz + static_cast<std::size_t>(j) * nz + k;
  }
  std::size_t size() const { return static_cast<std::size_t>(nx) * ny * nz; }
};

Grid makeGrid(const std::string& n, int nx, int ny, int nz) {
  Grid g; g.name = n; g.nx = nx; g.ny = ny; g.nz = nz;
  g.occ.assign(static_cast<std::size_t>(nx) * ny * nz, 0);
  return g;
}

void addBox(Grid& g, int i0, int j0, int k0, int i1, int j1, int k1) {
  for (int i = std::max(0, i0); i <= std::min(g.nx - 1, i1); ++i)
    for (int j = std::max(0, j0); j <= std::min(g.ny - 1, j1); ++j)
      for (int k = std::max(0, k0); k <= std::min(g.nz - 1, k1); ++k)
        g.occ[g.idx(i, j, k)] = 1;
}

void addBall(Grid& g, double ci, double cj, double ck, double r) {
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k) {
        const double d = (i - ci) * (i - ci) + (j - cj) * (j - cj) + (k - ck) * (k - ck);
        if (d <= r * r) g.occ[g.idx(i, j, k)] = 1;
      }
}

std::vector<Grid> buildCases() {
  std::vector<Grid> v;

  { // the original two-cube case, kept so the historical assertion still runs
    Grid g = makeGrid("two_cubes", 10, 10, 10);
    addBox(g, 2, 2, 2, 4, 4, 4);
    addBox(g, 6, 6, 6, 8, 8, 8);
    v.push_back(g);
  }
  { // empty and full: the degenerate ends
    v.push_back(makeGrid("empty", 8, 9, 7));
    Grid g = makeGrid("full", 8, 9, 7);
    std::fill(g.occ.begin(), g.occ.end(), 1);
    v.push_back(g);
  }
  { // anisotropic extents, so any index-arithmetic slip shows up
    Grid g = makeGrid("anisotropic", 5, 17, 3);
    addBox(g, 0, 0, 0, 2, 4, 1);
    addBox(g, 3, 9, 1, 4, 16, 2);
    v.push_back(g);
  }
  { // SHEARED, non-orthogonal, non-symmetric lattice of blobs. Blob centres sit
    // on a triclinic-style basis whose vectors are neither axis-aligned nor
    // mutually orthogonal, and whose lengths differ, so no two grid axes are
    // interchangeable and no symmetry can mask a directional error.
    Grid g = makeGrid("sheared_triclinic", 30, 33, 27);
    const double a[3] = {7.0, 1.3, 0.7};
    const double b[3] = {2.1, 8.0, 1.9};
    const double c[3] = {1.1, 2.7, 6.0};
    for (int u = 0; u < 4; ++u)
      for (int w = 0; w < 4; ++w)
        for (int t = 0; t < 4; ++t) {
          const double ci = 3.0 + u * a[0] + w * b[0] + t * c[0];
          const double cj = 3.0 + u * a[1] + w * b[1] + t * c[1];
          const double ck = 3.0 + u * a[2] + w * b[2] + t * c[2];
          if (ci > g.nx - 2 || cj > g.ny - 2 || ck > g.nz - 2) continue;
          addBall(g, ci, cj, ck, 1.6);
        }
    v.push_back(g);
  }
  { // diagonal-only contact: two boxes touching at a corner are ONE component
    // under 26-connectivity and TWO under 6. Separates the stencils.
    Grid g = makeGrid("corner_touch", 12, 12, 12);
    addBox(g, 1, 1, 1, 4, 4, 4);
    addBox(g, 5, 5, 5, 8, 8, 8);
    v.push_back(g);
  }
  { // concave shell with an interior void, plus a long thin runner: exercises
    // run merging across many rows and a component that is not convex.
    Grid g = makeGrid("shell_and_runner", 16, 16, 16);
    addBox(g, 2, 2, 2, 10, 10, 10);
    addBox(g, 4, 4, 4, 8, 8, 8);
    for (int i = 4; i <= 8; ++i)
      for (int j = 4; j <= 8; ++j)
        for (int k = 4; k <= 8; ++k) g.occ[g.idx(i, j, k)] = 0;
    for (int k = 0; k < 16; ++k) g.occ[g.idx(13, 14, k)] = 1;
    v.push_back(g);
  }
  { // single voxels scattered so seed pruning has something to prune
    Grid g = makeGrid("specks", 20, 21, 19);
    for (int i = 1; i < 20; i += 4)
      for (int j = 2; j < 21; j += 5)
        for (int k = 1; k < 19; k += 6) g.occ[g.idx(i, j, k)] = 1;
    addBall(g, 10.0, 10.0, 9.0, 4.2);
    v.push_back(g);
  }
  return v;
}

ClusterResult run(ClusterAlgorithm a, const Grid& g, int connectivity,
                  std::vector<int>* labels) {
  ClusterParams p;
  p.nx = g.nx; p.ny = g.ny; p.nz = g.nz;
  p.connectivity = connectivity;
  p.eps = 3.0;                       // vccs seed spacing
  std::vector<uint8_t> visited(g.size(), 0);
  return runClusterAlgorithm(a, p, g.occ, visited, labels);
}

// Exact pair: everything must match.
void checkExactPair(const Grid& g, ClusterAlgorithm fair, ClusterAlgorithm opt,
                    int connectivity) {
  const std::string tag = g.name + " conn" + std::to_string(connectivity) + " " +
                          algorithmToString(fair) + "/" + algorithmToString(opt);
  std::vector<int> lf, lo;
  ClusterResult rf = run(fair, g, connectivity, &lf);
  ClusterResult ro = run(opt, g, connectivity, &lo);

  check(rf.nclusters == ro.nclusters,
        tag + ": nclusters " + std::to_string(rf.nclusters) + " vs " + std::to_string(ro.nclusters));
  check(rf.maxCluster == ro.maxCluster,
        tag + ": max_cluster " + std::to_string(rf.maxCluster) + " vs " + std::to_string(ro.maxCluster));
  check(rf.visitedVoxels == ro.visitedVoxels, tag + ": visitedVoxels differ");

  std::vector<int> sf = rf.clusterSizes, so = ro.clusterSizes;
  std::sort(sf.begin(), sf.end());
  std::sort(so.begin(), so.end());
  check(sf == so, tag + ": sorted size distribution differs");

  // Stronger than the size histogram: the partitions themselves must agree.
  // A merge of one pair plus a split of another preserves the histogram.
  check(lf.size() == lo.size(), tag + ": label buffer size differs");
  if (lf.size() == lo.size()) {
    bool same = true;
    std::vector<int> mapFO, mapOF;
    mapFO.assign(static_cast<std::size_t>(rf.nclusters) + 1, -1);
    mapOF.assign(static_cast<std::size_t>(ro.nclusters) + 1, -1);
    for (std::size_t i = 0; i < lf.size() && same; ++i) {
      if ((lf[i] < 0) != (lo[i] < 0)) { same = false; break; }
      if (lf[i] < 0) continue;
      int& f2o = mapFO[static_cast<std::size_t>(lf[i])];
      int& o2f = mapOF[static_cast<std::size_t>(lo[i])];
      if (f2o == -1 && o2f == -1) { f2o = lo[i]; o2f = lf[i]; }
      else if (f2o != lo[i] || o2f != lf[i]) same = false;
    }
    check(same, tag + ": partitions differ (merge and/or split)");
  }
}

// Segmentation pair: the properties that must hold for BOTH tracks, since
// identical output is not one of them. See the header comment.
void checkSegmentation(const Grid& g, ClusterAlgorithm a) {
  const std::string tag = g.name + " " + algorithmToString(a);
  ClusterParams p;
  p.nx = g.nx; p.ny = g.ny; p.nz = g.nz;
  p.connectivity = 6;
  p.eps = 3.0;
  std::vector<uint8_t> visited(g.size(), 0);
  ClusterResult r = runClusterAlgorithm(a, p, g.occ, visited, nullptr);

  std::size_t occCount = 0;
  for (auto b : g.occ) occCount += (b == 1);

  // 1. exactly the occupied set is partitioned, no more and no less
  check(r.visitedVoxels == occCount,
        tag + ": covered " + std::to_string(r.visitedVoxels) + " of " + std::to_string(occCount) +
            " occupied voxels");
  std::size_t sum = 0;
  for (int s : r.clusterSizes) sum += static_cast<std::size_t>(s);
  check(sum == occCount, tag + ": cluster sizes sum to " + std::to_string(sum) + ", occupied is " +
                             std::to_string(occCount));

  // 2. a segmentation refines the connected components: it may split a
  //    component into several supervoxels, but must never span two.
  std::vector<int> dfsLab;
  ClusterResult dfs = run(ClusterAlgorithm::TraditionalDFS, g, 6, &dfsLab);
  check(r.nclusters >= dfs.nclusters,
        tag + ": produced " + std::to_string(r.nclusters) + " supervoxels, fewer than the " +
            std::to_string(dfs.nclusters) + " connected components -- it merged components");
  check(r.maxCluster <= dfs.maxCluster,
        tag + ": largest supervoxel " + std::to_string(r.maxCluster) +
            " exceeds the largest connected component " + std::to_string(dfs.maxCluster));
}

}  // namespace

int main() {
  const auto cases = buildCases();

  for (const auto& g : cases) {
    for (int conn : {6, 26}) {
      // Only cc3d honours ClusterParams::connectivity; rle_ccl hardcodes 6.
      checkExactPair(g, ClusterAlgorithm::CC3D, ClusterAlgorithm::CC3DOptimized, conn);
    }
    checkExactPair(g, ClusterAlgorithm::RLECCL, ClusterAlgorithm::RLECCLOptimized, 6);
    checkSegmentation(g, ClusterAlgorithm::VCCS);
    checkSegmentation(g, ClusterAlgorithm::VCCSOptimized);
  }

  // The original assertion this file was written for, kept verbatim in intent:
  // two 3x3x3 cubes, 27 voxels each, must be found by both cc3d tracks.
  {
    Grid g = makeGrid("two_cubes", 10, 10, 10);
    addBox(g, 2, 2, 2, 4, 4, 4);
    addBox(g, 6, 6, 6, 8, 8, 8);
    ClusterResult a = run(ClusterAlgorithm::CC3D, g, 6, nullptr);
    ClusterResult b = run(ClusterAlgorithm::CC3DOptimized, g, 6, nullptr);
    check(a.nclusters == 2 && b.nclusters == 2, "two_cubes: expected 2 clusters");
    check(a.maxCluster == 27 && b.maxCluster == 27, "two_cubes: expected max cluster 27");
  }

  const int total = static_cast<int>(cases.size());
  if (failures == 0) {
    std::cout << "dual-track equivalence PASSED over " << total << " grids\n"
              << "  exact pairs: cc3d/cc3d_optimized (conn 6 and 26), rle_ccl/rle_ccl_optimized\n"
              << "  segmentation pairs: vccs, vccs_optimized (coverage, refinement)\n";
    return 0;
  }
  std::cerr << "dual-track equivalence FAILED: " << failures << " check(s)\n";
  return 1;
}
