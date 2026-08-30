// Task 4: differential testing of the algorithms that claim exact labelling.
//
// TraditionalDFS, CC3D, CC3DOptimized, GCBD, RLE-CCL and (since Task 12)
// RLE-CCL-optimized all claim to compute connected components exactly. If that
// is true they must agree with each other on every input, and all of them must
// agree with a straightforward reference. Octree-CCL sat in this list until
// Task 12 removed it. RLE-CCL and its optimized variant are structurally clever
// (run-merging; runs rather than voxels as the union-find domain) in ways that
// make a merge-plus-split error easy to introduce and invisible to a count-only
// comparison.
//
// So this compares three things of increasing strength:
//   1. component count
//   2. sorted size histogram
//   3. the full label partition, up to relabelling
// A method that merges one pair of components and splits another passes (1),
// can pass (2), and cannot pass (3). Only (3) actually tests the claim.
//
// The verdict is taken against an INDEPENDENT reference BFS written here, not
// only by majority vote: six implementations that share a misconception would
// out-vote a correct one. Majority vote is still computed and reported, since
// it is what identifies the outlier when the reference itself is in doubt.
//
// The PBC on/off axis originally planned for this task does not exist -- no
// algorithm in this codebase has any periodic boundary handling -- and is
// replaced by a connectivity sweep. Note that only cc3d and cc3dOptimized read
// ClusterParams::connectivity; gcbd and both rleCCL variants hardcode
// 6-connectivity and traditionalDFS hardcodes deltas6. At connectivity 26 the
// group is therefore NOT computing the same thing, and diffing all of them
// would report spurious failures. Only the two that honour the parameter are
// compared there, against a 26-connected reference.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "cluster/Algorithms.hpp"

using bls::ClusterAlgorithm;
using bls::ClusterParams;
using bls::ClusterResult;

namespace {

struct Grid3 {
  int nx{0}, ny{0}, nz{0};
  std::vector<uint8_t> occ;

  std::size_t index(int i, int j, int k) const {
    return static_cast<std::size_t>(i) * ny * nz + static_cast<std::size_t>(j) * nz + k;
  }
  std::size_t size() const { return static_cast<std::size_t>(nx) * ny * nz; }
  std::size_t occupied() const {
    std::size_t n = 0;
    for (uint8_t v : occ)
      if (v) ++n;
    return n;
  }
};

Grid3 makeGrid(int nx, int ny, int nz) {
  Grid3 g;
  g.nx = nx; g.ny = ny; g.nz = nz;
  g.occ.assign(g.size(), 0);
  return g;
}

const int kD6[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

// Independent reference: plain BFS over the occupancy grid. Deliberately the
// dumbest correct implementation -- no union-find, no runs, no decision tree -- so it
// shares no machinery, and therefore no potential misconception, with any
// method under test.
std::vector<int> referenceLabels(const Grid3& g, int connectivity) {
  std::vector<std::array<int, 3>> dirs;
  if (connectivity == 26) {
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
          if (dx || dy || dz) dirs.push_back({dx, dy, dz});
  } else {
    for (const auto& d : kD6) dirs.push_back({d[0], d[1], d[2]});
  }

  std::vector<int> labels(g.size(), -1);
  std::vector<std::size_t> stack;
  int next = 0;
  for (int i = 0; i < g.nx; ++i) {
    for (int j = 0; j < g.ny; ++j) {
      for (int k = 0; k < g.nz; ++k) {
        std::size_t s = g.index(i, j, k);
        if (!g.occ[s] || labels[s] >= 0) continue;
        const int id = next++;
        stack.clear();
        stack.push_back(s);
        labels[s] = id;
        while (!stack.empty()) {
          std::size_t cur = stack.back();
          stack.pop_back();
          int ci = static_cast<int>(cur / (static_cast<std::size_t>(g.ny) * g.nz));
          int rem = static_cast<int>(cur % (static_cast<std::size_t>(g.ny) * g.nz));
          int cj = rem / g.nz, ck = rem % g.nz;
          for (const auto& d : dirs) {
            int ni = ci + d[0], nj = cj + d[1], nk = ck + d[2];
            if (ni < 0 || nj < 0 || nk < 0 || ni >= g.nx || nj >= g.ny || nk >= g.nz) continue;
            std::size_t nidx = g.index(ni, nj, nk);
            if (!g.occ[nidx] || labels[nidx] >= 0) continue;
            labels[nidx] = id;
            stack.push_back(nidx);
          }
        }
      }
    }
  }
  return labels;
}

// Two label arrays describe the same partition iff they are equal after
// renumbering in order of first appearance on a linear scan. This is the
// harness's job precisely so no algorithm has to pay for a canonical ordering
// on its timed path.
std::vector<int> canonicalize(const std::vector<int>& labels) {
  std::vector<int> out(labels.size(), -1);
  std::map<int, int> remap;
  int next = 0;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    if (labels[i] < 0) continue;
    auto it = remap.find(labels[i]);
    if (it == remap.end()) it = remap.emplace(labels[i], next++).first;
    out[i] = it->second;
  }
  return out;
}

std::vector<int> sizeHistogram(const std::vector<int>& canonical) {
  std::map<int, int> counts;
  for (int v : canonical)
    if (v >= 0) ++counts[v];
  std::vector<int> h;
  for (const auto& kv : counts) h.push_back(kv.second);
  std::sort(h.begin(), h.end(), std::greater<int>());
  return h;
}

int componentCount(const std::vector<int>& canonical) {
  int mx = -1;
  for (int v : canonical) mx = std::max(mx, v);
  return mx + 1;
}

// --- grid generators --------------------------------------------------------

std::vector<std::pair<std::string, Grid3>> structuredCases() {
  std::vector<std::pair<std::string, Grid3>> out;
  auto add = [&](const std::string& n, Grid3 g) { out.emplace_back(n, std::move(g)); };

  for (int n : {4, 5, 8, 13}) {
    { Grid3 g = makeGrid(n, n, n); add("empty-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n); g.occ[g.index(n/2, n/2, n/2)] = 1;
      add("single-voxel-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n); std::fill(g.occ.begin(), g.occ.end(), 1);
      add("full-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n);  // a plane through the middle
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) g.occ[g.index(i, j, n/2)] = 1;
      add("plane-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n);  // a line along x
      for (int i = 0; i < n; ++i) g.occ[g.index(i, n/2, n/2)] = 1;
      add("line-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n);  // body diagonal: disconnected under 6, one under 26
      for (int i = 0; i < n; ++i) g.occ[g.index(i, i, i)] = 1;
      add("diagonal-chain-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n);  // two blocks meeting at exactly one corner
      const int h = n / 2;
      for (int i = 0; i < h; ++i) for (int j = 0; j < h; ++j) for (int k = 0; k < h; ++k)
        g.occ[g.index(i, j, k)] = 1;
      for (int i = h; i < n; ++i) for (int j = h; j < n; ++j) for (int k = h; k < n; ++k)
        g.occ[g.index(i, j, k)] = 1;
      add("corner-touching-" + std::to_string(n), g); }
    { Grid3 g = makeGrid(n, n, n);  // spanning cluster: a full face-to-face column
      for (int k = 0; k < n; ++k) { g.occ[g.index(0, 0, k)] = 1; g.occ[g.index(n-1, n-1, k)] = 1; }
      add("spanning-" + std::to_string(n), g); }
  }

  // Nested shells: concentric hollow cubes, each a separate component, which
  // catches an algorithm that leaks across an enclosed cavity.
  for (int n : {9, 13, 17}) {
    Grid3 g = makeGrid(n, n, n);
    for (int r = 0; r <= n / 2; r += 2) {
      for (int i = r; i < n - r; ++i)
        for (int j = r; j < n - r; ++j)
          for (int k = r; k < n - r; ++k) {
            const bool onShell = i == r || j == r || k == r || i == n-1-r || j == n-1-r || k == n-1-r;
            if (onShell) g.occ[g.index(i, j, k)] = 1;
          }
    }
    g = g;
    out.emplace_back("nested-shells-" + std::to_string(n), g);
  }

  // Two clusters separated by a single empty plane, in each orientation.
  for (int axis = 0; axis < 3; ++axis) {
    const int n = 11;
    Grid3 g = makeGrid(n, n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k) {
          const int c = axis == 0 ? i : (axis == 1 ? j : k);
          if (c != n / 2) g.occ[g.index(i, j, k)] = 1;
        }
    out.emplace_back("split-plane-axis" + std::to_string(axis), g);
  }

  // Non-cubic shapes, so no algorithm passes by assuming nx==ny==nz.
  {
    Grid3 g = makeGrid(3, 17, 29);
    std::mt19937_64 rng(9);
    std::uniform_real_distribution<double> u(0, 1);
    for (auto& v : g.occ) v = u(rng) < 0.31 ? 1 : 0;
    out.emplace_back("slab-3x17x29", g);
  }
  {
    Grid3 g = makeGrid(29, 3, 17);
    std::mt19937_64 rng(10);
    std::uniform_real_distribution<double> u(0, 1);
    for (auto& v : g.occ) v = u(rng) < 0.31 ? 1 : 0;
    out.emplace_back("slab-29x3x17", g);
  }
  {
    Grid3 g = makeGrid(1, 1, 64);
    for (int k = 0; k < 64; k += 2) g.occ[g.index(0, 0, k)] = 1;
    out.emplace_back("comb-1x1x64", g);
  }
  return out;
}

// --- minimisation -----------------------------------------------------------

// Greedy delta-debugging: clear occupied voxels one at a time, keeping a
// removal only if the disagreement survives it, then crop to the bounding box
// of what remains. Reports the smallest grid that still reproduces.
template <typename FailsFn>
Grid3 shrink(Grid3 g, const FailsFn& fails) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t i = 0; i < g.occ.size(); ++i) {
      if (!g.occ[i]) continue;
      g.occ[i] = 0;
      if (fails(g)) {
        changed = true;
      } else {
        g.occ[i] = 1;
      }
    }
  }
  int lo[3] = {g.nx, g.ny, g.nz}, hi[3] = {-1, -1, -1};
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k)
        if (g.occ[g.index(i, j, k)]) {
          int c[3] = {i, j, k};
          for (int a = 0; a < 3; ++a) { lo[a] = std::min(lo[a], c[a]); hi[a] = std::max(hi[a], c[a]); }
        }
  if (hi[0] < 0) return g;
  Grid3 cropped = makeGrid(hi[0]-lo[0]+1, hi[1]-lo[1]+1, hi[2]-lo[2]+1);
  for (int i = lo[0]; i <= hi[0]; ++i)
    for (int j = lo[1]; j <= hi[1]; ++j)
      for (int k = lo[2]; k <= hi[2]; ++k)
        if (g.occ[g.index(i, j, k)])
          cropped.occ[cropped.index(i-lo[0], j-lo[1], k-lo[2])] = 1;
  return fails(cropped) ? cropped : g;
}

void printGrid(const Grid3& g) {
  std::printf("      %dx%dx%d, %zu occupied\n", g.nx, g.ny, g.nz, g.occupied());
  if (g.size() > 4096) { std::printf("      (too large to print)\n"); return; }
  for (int i = 0; i < g.nx; ++i) {
    std::printf("      i=%d\n", i);
    for (int j = 0; j < g.ny; ++j) {
      std::printf("        ");
      for (int k = 0; k < g.nz; ++k) std::printf("%c", g.occ[g.index(i, j, k)] ? '#' : '.');
      std::printf("\n");
    }
  }
}

// --- the comparison itself --------------------------------------------------

struct Method {
  ClusterAlgorithm algo;
  const char* name;
};

const Method kSix[] = {
    {ClusterAlgorithm::TraditionalDFS, "traditional_dfs"},
    {ClusterAlgorithm::CC3D,           "cc3d"},
    {ClusterAlgorithm::CC3DOptimized,  "cc3d_optimized"},
    {ClusterAlgorithm::GCBD,           "gcbd"},
    {ClusterAlgorithm::RLECCL,         "rle_ccl"},
    {ClusterAlgorithm::RLECCLOptimized,"rle_ccl_optimized"},
};

struct Tally {
  long countMismatch = 0;
  long histMismatch = 0;
  long partitionMismatch = 0;
  bool reported = false;
};

std::vector<int> runOne(const Method& m, const Grid3& g, int connectivity, ClusterResult* out) {
  ClusterParams p;
  p.nx = g.nx; p.ny = g.ny; p.nz = g.nz;
  p.connectivity = connectivity;
  std::vector<uint8_t> visited(g.size(), 0);
  std::vector<int> labels;
  ClusterResult r = bls::runClusterAlgorithm(m.algo, p, g.occ, visited, &labels);
  if (out) *out = r;
  return canonicalize(labels);
}

}  // namespace

// Positive control for the comparison machinery itself.
//
// A harness that compares nothing reports zero mismatches, which looks exactly
// like success. So before trusting a clean run, inject the two faults this task
// exists to catch and confirm each is detected -- including the merge-plus-split
// that leaves the component COUNT correct, which is the whole reason the
// partition comparison is here.
bool selfCheck() {
  Grid3 g = makeGrid(7, 7, 7);
  // Four separated 2x2x2 blocks, so there is something to merge and split.
  const int origins[4][3] = {{0,0,0},{4,0,0},{0,4,0},{4,4,4}};
  for (const auto& o : origins)
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
      g.occ[g.index(o[0]+a, o[1]+b, o[2]+c)] = 1;

  const std::vector<int> ref = canonicalize(referenceLabels(g, 6));
  const int refCount = componentCount(ref);
  const std::vector<int> refHist = sizeHistogram(ref);
  if (refCount != 4) {
    std::printf("  SELF-CHECK FAILED: reference found %d components, expected 4\n", refCount);
    return false;
  }

  bool ok = true;

  // Fault A: merge components 0 and 1. Count drops, so even a count-only
  // comparison sees this one.
  {
    std::vector<int> bad = ref;
    for (int& v : bad) if (v == 1) v = 0;
    const std::vector<int> c = canonicalize(bad);
    if (componentCount(c) == refCount || c == ref) {
      std::printf("  SELF-CHECK FAILED: a merged partition was not detected\n");
      ok = false;
    }
  }

  // Fault B: merge 0 with 1 AND split 2 in half. Component count is back to 4
  // and the size histogram is {16,4,4,8} sorted -> {16,8,4,4} vs reference
  // {8,8,8,8}: the count check passes, the histogram catches it here, and the
  // partition check catches it unconditionally. This is the defect shape that a
  // count-only diff would have missed entirely.
  {
    std::vector<int> bad = ref;
    for (int& v : bad) if (v == 1) v = 0;
    int seen = 0;
    for (int& v : bad) if (v == 2 && seen++ >= 4) v = 99;
    const std::vector<int> c = canonicalize(bad);
    if (componentCount(c) != refCount) {
      std::printf("  SELF-CHECK FAILED: merge+split changed the count, so it is not the "
                  "count-preserving case this check is for\n");
      ok = false;
    }
    if (c == ref) {
      std::printf("  SELF-CHECK FAILED: a merge+split partition compared equal to the reference\n");
      ok = false;
    }
    if (sizeHistogram(c) == refHist) {
      std::printf("  SELF-CHECK FAILED: merge+split produced an identical size histogram\n");
      ok = false;
    }
  }

  // Fault C: relabelling alone must NOT be reported as a difference, or every
  // comparison below would be a false positive.
  {
    std::vector<int> permuted = ref;
    for (int& v : permuted) if (v >= 0) v = 1000 - v;
    if (canonicalize(permuted) != ref) {
      std::printf("  SELF-CHECK FAILED: pure relabelling was reported as a partition change\n");
      ok = false;
    }
  }

  // Fault D: the minimiser must actually shrink, and its result must still fail.
  {
    auto fails = [](const Grid3& cand) { return cand.occupied() >= 2; };
    Grid3 small = shrink(g, fails);
    if (!fails(small) || small.occupied() >= g.occupied()) {
      std::printf("  SELF-CHECK FAILED: minimiser did not shrink (%zu -> %zu occupied)\n",
                  g.occupied(), small.occupied());
      ok = false;
    }
  }

  std::printf("  self-check: %s\n", ok ? "comparison detects merge, split, merge+split, and "
                                         "minimises; pure relabelling is not flagged"
                                       : "FAILED");
  return ok;
}

int main(int argc, char** argv) {
  const long targetRandom = argc > 1 ? std::stol(argv[1]) : 5000;

  std::printf("== harness self-check ==\n");
  if (!selfCheck()) {
    std::printf("\nAborting: the comparison cannot be trusted, so a clean run would be "
                "meaningless.\n");
    return 2;
  }
  std::printf("\n");

  std::map<std::string, Tally> tally;
  for (const auto& m : kSix) tally[m.name];

  long gridsChecked = 0, comparisons = 0, disagreeingGrids = 0;
  std::size_t voxelsChecked = 0;

  // Checks one grid at one connectivity. Returns true if anything disagreed.
  auto checkGrid = [&](const Grid3& g, const std::string& name, int connectivity) {
    const std::vector<int> refCanon = canonicalize(referenceLabels(g, connectivity));
    const int refCount = componentCount(refCanon);
    const std::vector<int> refHist = sizeHistogram(refCanon);

    // At connectivity 26 only the two methods that actually read the parameter
    // are comparable; see the header comment.
    const bool only26Aware = connectivity == 26;

    std::vector<std::string> partitionOutliers;
    std::map<std::vector<int>, std::vector<std::string>> byPartition;
    bool any = false;

    for (const auto& m : kSix) {
      if (only26Aware && m.algo != ClusterAlgorithm::CC3D &&
          m.algo != ClusterAlgorithm::CC3DOptimized) {
        continue;
      }
      ClusterResult r;
      const std::vector<int> canon = runOne(m, g, connectivity, &r);
      ++comparisons;
      byPartition[canon].push_back(m.name);

      auto& t = tally[m.name];
      bool bad = false;
      if (r.nclusters != refCount) { ++t.countMismatch; bad = true; }
      std::vector<int> reported = r.clusterSizes;
      std::sort(reported.begin(), reported.end(), std::greater<int>());
      if (reported != refHist) { ++t.histMismatch; bad = true; }
      if (canon != refCanon) { ++t.partitionMismatch; bad = true; partitionOutliers.push_back(m.name); }

      if (bad) {
        any = true;
        if (!t.reported) {
          t.reported = true;
          std::printf("\n  DISAGREEMENT  %s  on '%s'  connectivity=%d\n", m.name, name.c_str(),
                      connectivity);
          std::printf("      reference: %d components, %zu occupied\n", refCount, g.occupied());
          std::printf("      %-16s %d components\n", m.name, r.nclusters);

          // Minimise, using this method's disagreement with the reference as
          // the predicate.
          Method mm = m;
          const int conn = connectivity;
          Grid3 small = shrink(g, [&](const Grid3& cand) {
            if (cand.occupied() == 0) return false;
            const std::vector<int> rc = canonicalize(referenceLabels(cand, conn));
            ClusterResult rr;
            const std::vector<int> cc = runOne(mm, cand, conn, &rr);
            return cc != rc || rr.nclusters != componentCount(rc);
          });
          std::printf("      minimal reproducing grid:\n");
          printGrid(small);
        }
      }
    }

    if (any) {
      ++disagreeingGrids;
      // Majority vote, as specified: whichever partition the fewest methods
      // produce is the outlier. Reported alongside the reference verdict, which
      // is the one actually trusted.
      if (byPartition.size() > 1) {
        std::size_t best = 0;
        for (const auto& kv : byPartition) best = std::max(best, kv.second.size());
        std::printf("      majority vote: ");
        for (const auto& kv : byPartition) {
          if (kv.second.size() < best) {
            for (const auto& n : kv.second) std::printf("%s ", n.c_str());
          }
        }
        std::printf("in the minority\n");
      } else {
        std::printf("      majority vote: unanimous among the methods; "
                    "they disagree with the independent reference\n");
      }
    }
    ++gridsChecked;
    voxelsChecked += g.size();
    return any;
  };

  // ---- structured cases ----
  std::printf("== structured cases ==\n");
  const auto structured = structuredCases();
  for (const auto& sc : structured) {
    checkGrid(sc.second, sc.first, 6);
    checkGrid(sc.second, sc.first, 26);
  }
  std::printf("  %zu structured grids x 2 connectivities\n", structured.size());

  // ---- random grids ----
  // Sizes 4^3..64^3, occupancy 0.01..0.95 with the site-percolation threshold
  // for the simple cubic lattice (p_c ~= 0.3116) sampled densely, since that is
  // where the component structure is most fragile and a merge or split is most
  // likely to be reachable.
  std::printf("\n== random grids ==\n");
  const int sizes[] = {4, 6, 8, 12, 16, 24, 32, 48, 64};
  const int weights[] = {8, 8, 12, 14, 16, 16, 14, 8, 4};
  const double occs[] = {0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25,
                         0.28, 0.30, 0.3116, 0.32, 0.34, 0.40, 0.50, 0.65, 0.80, 0.95};

  std::mt19937_64 rng(20260829ULL);
  std::discrete_distribution<int> sizePick(std::begin(weights), std::end(weights));
  std::uniform_int_distribution<int> occPick(0, static_cast<int>(std::size(occs)) - 1);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::uniform_int_distribution<int> connPick(0, 9);

  for (long t = 0; t < targetRandom; ++t) {
    const int n = sizes[sizePick(rng)];
    // Occasionally make it non-cubic, so shape assumptions are exercised too.
    int nx = n, ny = n, nz = n;
    if (u01(rng) < 0.2) {
      std::uniform_int_distribution<int> jitter(1, std::max(1, n));
      nx = jitter(rng); ny = jitter(rng); nz = jitter(rng);
    }
    const double p = occs[occPick(rng)];
    Grid3 g = makeGrid(nx, ny, nz);
    for (auto& v : g.occ) v = u01(rng) < p ? 1 : 0;

    // One in ten grids is also checked at connectivity 26 (two methods only).
    const int conn = connPick(rng) == 0 ? 26 : 6;
    checkGrid(g, "random-" + std::to_string(nx) + "x" + std::to_string(ny) + "x" +
                     std::to_string(nz) + "-p" + std::to_string(p) + "-#" + std::to_string(t),
              conn);
  }

  // ---- summary ----
  std::printf("\n== summary ==\n");
  std::printf("  grids checked      : %ld\n", gridsChecked);
  std::printf("  voxels checked     : %zu\n", voxelsChecked);
  std::printf("  algorithm runs     : %ld\n", comparisons);
  std::printf("  grids with any disagreement: %ld\n", disagreeingGrids);
  std::printf("\n  %-18s %12s %12s %12s\n", "algorithm", "count", "histogram", "partition");
  long total = 0;
  for (const auto& m : kSix) {
    const auto& t = tally[m.name];
    std::printf("  %-18s %12ld %12ld %12ld\n", m.name, t.countMismatch, t.histMismatch,
                t.partitionMismatch);
    total += t.countMismatch + t.histMismatch + t.partitionMismatch;
  }
  std::printf("\n  total mismatches: %ld\n", total);
  return total == 0 ? 0 : 1;
}
