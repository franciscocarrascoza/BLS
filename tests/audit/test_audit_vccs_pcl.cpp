// Created: 2026-09-26T18:10+02:00 | by: CC audit (manager) | purpose: WP4 tests — VCCS in its published (PCL) form (decision D6; findings VCCS-1/2/5/6)
//
// Guards finding / decision IDs:
//   VCCS-REF   vccs_optimized == an independent, deliberately plain re-implementation of PCL's
//              SupervoxelClustering on the grid (full-grid arrays, no binary search, no sets): same
//              cluster sizes, seed candidates, seeds kept and threshold, on random grids, S = 2..6,
//              connectivity 6 and 26, PBC none and PBC xyz
//   VCCS-SEED  every candidate is the occupied voxel nearest its seed-cell centre (brute force over
//              all occupied voxels), one per occupied cell, duplicates merged
//   VCCS-PRUNE threshold = 0.05*pi*(S/2)^2 (PCL), count = #{|v-c|^2 < (S/2)^2} (FLANN strict), kept iff >
//   VCCS-26    growth and the remainder pass use the method's connectivity: a diagonal-only chain is one
//              cluster at 26 (PCL) and seven at 6
//   VCCS-TRANS PBC xyz with every N a multiple of S: rolling the grid by multiples of S leaves the seed-cell
//              count unchanged. The cluster sizes are NOT roll-invariant, and are not expected to be:
//              PCL's growth is order-dependent (supervoxels expand one after another within a round and
//              a later one can take a voxel from an earlier one), and the order follows the scan. The
//              size of that effect is printed as a measurement.
//   VCCS-COVER every occupied voxel is counted exactly once (sum of sizes = occupied voxels)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "audit_common.hpp"

namespace {

// Plain reference of PCL's algorithm (see Algorithms.cpp vccsOptimized for the rules).
struct RefResult {
  std::vector<int> sizes;  // descending
  int candidates{0}, kept{0};
  double threshold{0};
};

RefResult pclReference(const audit::Grid3& g, int S, bool periodic, int conn) {
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  auto wrap = [](int v, int n) { return ((v % n) + n) % n; };
  auto occAt = [&](int i, int j, int k, int& out) {
    if (periodic) { i = wrap(i, nx); j = wrap(j, ny); k = wrap(k, nz); }
    else if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return false;
    out = static_cast<int>(g.index(i, j, k));
    return g.occ[static_cast<std::size_t>(out)] != 0;
  };
  auto d1 = [&](double a, double b, int n) {
    double d = a - b;
    if (periodic) d -= n * std::round(d / n);
    return d;
  };
  auto coords = [&](int v, int& i, int& j, int& k) { i = v / (ny * nz); j = (v / nz) % ny; k = v % nz; };
  auto d2 = [&](int v, double x, double y, double z) {
    int i, j, k; coords(v, i, j, k);
    const double a = d1(i, x, nx), b = d1(j, y, ny), c = d1(k, z, nz);
    return a * a + b * b + c * c;
  };
  std::vector<int> occ;
  for (std::size_t v = 0; v < g.size(); ++v) if (g.occ[v]) occ.push_back(static_cast<int>(v));
  RefResult r;
  // seeding: brute force nearest to each occupied cell centre
  std::map<std::array<int, 3>, int> cellBest;
  const double h = 0.5 * (S - 1);
  for (int v : occ) {
    int i, j, k; coords(v, i, j, k);
    cellBest.emplace(std::array<int, 3>{i / S, j / S, k / S}, -1);
  }
  std::vector<int> cand;
  for (auto& kv : cellBest) {
    const double cx = kv.first[0] * S + h, cy = kv.first[1] * S + h, cz = kv.first[2] * S + h;
    int best = -1; double bd = 0;
    for (int v : occ) {
      const double d = d2(v, cx, cy, cz);
      if (best < 0 || d < bd) { best = v; bd = d; }  // occ ascending -> ties keep the smaller index
    }
    cand.push_back(best);
  }
  std::sort(cand.begin(), cand.end());
  cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
  r.candidates = static_cast<int>(cand.size());
  // pruning
  const double R = 0.5 * S;
  r.threshold = 0.05 * R * R * 3.1415926536;
  std::vector<int> seeds;
  for (int c : cand) {
    int i, j, k; coords(c, i, j, k);
    int num = 0;
    for (int v : occ) if (d2(v, i, j, k) < R * R) ++num;
    if (num > r.threshold) seeds.push_back(c);
  }
  r.kept = static_cast<int>(seeds.size());
  // growth
  const int n = static_cast<int>(seeds.size());
  std::vector<int> owner(g.size(), -1);
  std::vector<double> dist(g.size(), std::numeric_limits<double>::max());
  std::vector<std::vector<int>> leaves(n);  // kept sorted
  std::vector<std::array<double, 3>> cen(n);
  std::vector<bool> alive(n, true);
  for (int s = 0; s < n; ++s) {
    owner[seeds[s]] = s; leaves[s] = {seeds[s]};
    int i, j, k; coords(seeds[s], i, j, k); cen[s] = {double(i), double(j), double(k)};
  }
  const int depth = static_cast<int>(1.8 * S);
  for (int round = 1; round < depth; ++round) {
    for (int s = 0; s < n; ++s) {
      if (!alive[s]) continue;
      std::vector<int> fresh;
      for (int leaf : leaves[s]) {
        int i, j, k; coords(leaf, i, j, k);
        for (int a = -1; a <= 1; ++a) for (int b = -1; b <= 1; ++b) for (int c = -1; c <= 1; ++c) {
          const int m = std::abs(a) + std::abs(b) + std::abs(c);
          if (m == 0 || (conn == 6 && m != 1)) continue;
          int v;
          if (!occAt(i + a, j + b, k + c, v) || owner[v] == s) continue;
          const double d = std::sqrt(d2(v, cen[s][0], cen[s][1], cen[s][2]));
          if (d < dist[v]) {
            dist[v] = d;
            if (owner[v] >= 0) {
              auto& L = leaves[owner[v]];
              L.erase(std::find(L.begin(), L.end(), v));
            }
            owner[v] = s; fresh.push_back(v);
          }
        }
      }
      for (int v : fresh) leaves[s].insert(std::lower_bound(leaves[s].begin(), leaves[s].end(), v), v);
    }
    for (int s = 0; s < n; ++s) {
      if (!alive[s]) continue;
      if (leaves[s].empty()) { alive[s] = false; continue; }
      std::array<double, 3> c0 = cen[s], sum{0, 0, 0};
      for (int leaf : leaves[s]) {
        int i, j, k; coords(leaf, i, j, k);
        sum[0] += c0[0] + d1(i, c0[0], nx); sum[1] += c0[1] + d1(j, c0[1], ny); sum[2] += c0[2] + d1(k, c0[2], nz);
      }
      for (int a = 0; a < 3; ++a) cen[s][a] = sum[a] / leaves[s].size();
    }
  }
  std::vector<int> label(g.size(), -1);
  int next = 0;
  for (int s = 0; s < n; ++s) {
    if (!alive[s] || leaves[s].empty()) continue;
    for (int v : leaves[s]) label[v] = next;
    ++next;
  }
  for (int v0 : occ) {
    if (label[v0] >= 0) continue;
    std::vector<int> st{v0}; label[v0] = next;
    while (!st.empty()) {
      const int cur = st.back(); st.pop_back();
      int i, j, k; coords(cur, i, j, k);
      for (int a = -1; a <= 1; ++a) for (int b = -1; b <= 1; ++b) for (int c = -1; c <= 1; ++c) {
        const int m = std::abs(a) + std::abs(b) + std::abs(c);
        int v;
        if (m && (conn == 26 || m == 1) && occAt(i + a, j + b, k + c, v) && label[v] < 0) { label[v] = next; st.push_back(v); }
      }
    }
    ++next;
  }
  r.sizes = audit::sizesOf(label);
  return r;
}

bls::ClusterResult runVccs(const audit::Grid3& g, int S, bool periodic, int conn = 26) {
  auto p = audit::paramsFor(g);
  p.eps = S;
  p.connectivity = conn;
  p.periodic = periodic;
  return audit::runUnlabelled(bls::ClusterAlgorithm::VCCSOptimized, g, p);
}

void reference() {
  std::mt19937_64 rng(1301);
  std::uniform_int_distribution<int> dim(4, 18);
  int cases = 0;
  for (int S : {2, 3, 4, 5, 6})
    for (double p : {0.05, 0.2, 0.45, 0.8})
      for (bool periodic : {false, true})
        for (int conn : {6, 26})
        for (int rep = 0; rep < 4; ++rep) {
          const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
          const auto ref = pclReference(g, S, periodic, conn);
          const auto got = runVccs(g, S, periodic, conn);
          ++cases;
          const std::string tag = "S=" + std::to_string(S) + " conn=" + std::to_string(conn) + " p=" + std::to_string(p) + (periodic ? " PBC" : "") +
                                  " " + std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" + std::to_string(g.nz);
          AUDIT_CHECK("VCCS-REF", got.clusterSizes == ref.sizes,
                      tag + ": sizes differ (" + std::to_string(got.nclusters) + " vs " + std::to_string(ref.sizes.size()) + " clusters)");
          AUDIT_CHECK("VCCS-SEED", got.seedCandidates == ref.candidates,
                      tag + ": candidates " + std::to_string(got.seedCandidates) + " vs " + std::to_string(ref.candidates));
          AUDIT_CHECK("VCCS-PRUNE", got.seedsPlaced == ref.kept && std::fabs(got.seedPruneThreshold - ref.threshold) < 1e-9,
                      tag + ": kept " + std::to_string(got.seedsPlaced) + " vs " + std::to_string(ref.kept));
          long long total = 0;
          for (int sz : got.clusterSizes) total += sz;
          long long occ = 0;
          for (auto v : g.occ) occ += v;
          AUDIT_CHECK("VCCS-COVER", total == occ && got.visitedVoxels == static_cast<std::size_t>(occ), tag + ": coverage");
        }
  std::printf("measure: VCCS-REF %d random cases (S 2-6, p 0.05-0.8, connectivity 6/26, PBC none/xyz)\n", cases);
}

void adjacency26() {
  // A diagonal-only chain: 6-connected it is 7 components, 26-connected one.
  auto g = audit::emptyGrid(12, 12, 12);
  for (int t = 2; t < 9; ++t) g.set(t, t, t);
  const auto r26 = runVccs(g, 12, false, 26);  // S = 12: one seed cell
  const auto r6 = runVccs(g, 12, false, 6);
  AUDIT_CHECK("VCCS-26", r26.nclusters == 1, "26-adjacency: diagonal chain gave " + std::to_string(r26.nclusters) + " clusters");
  AUDIT_CHECK("VCCS-26", r6.nclusters == 7, "6-adjacency: diagonal chain gave " + std::to_string(r6.nclusters) + " clusters (want 7)");
}

void periodicRoll() {
  std::mt19937_64 rng(77);
  int changed = 0, total = 0, maxDelta = 0;
  for (int S : {2, 3, 4}) {
    for (int rep = 0; rep < 10; ++rep) {
      const auto g = audit::randomGrid(rng, 3 * S, 4 * S, 2 * S + (S == 2 ? 2 : 0), 0.3);
      auto roll = audit::emptyGrid(g.nx, g.ny, g.nz);
      const int ti = S, tj = 2 * S, tk = S;
      for (int i = 0; i < g.nx; ++i) for (int j = 0; j < g.ny; ++j) for (int k = 0; k < g.nz; ++k)
        if (g.at(i, j, k)) roll.set((i + ti) % g.nx, (j + tj) % g.ny, (k + tk) % g.nz);
      const auto a = runVccs(g, S, true), b = runVccs(roll, S, true);
      AUDIT_CHECK("VCCS-TRANS", a.seedCandidates == b.seedCandidates,
                  "S=" + std::to_string(S) + ": seed-cell count changes under a roll by multiples of S (" +
                      std::to_string(a.seedCandidates) + " vs " + std::to_string(b.seedCandidates) + ")");
      ++total;
      if (a.clusterSizes != b.clusterSizes) ++changed;
      maxDelta = std::max(maxDelta, std::abs(a.nclusters - b.nclusters));
    }
  }
  std::printf("measure: VCCS order dependence: %d of %d periodic rolls change the size multiset; "
              "largest change in cluster count %d\n", changed, total, maxDelta);
}

}  // namespace

int main() {
  reference();
  adjacency26();
  periodicRoll();
  return audit::finish("test_audit_vccs_pcl");
}
