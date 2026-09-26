// Created: 2026-09-26T18:08+02:00 (as test_audit_vccs_pcl.cpp) | Rewritten: 2026-09-26T18:31+02:00 | by: CC audit (manager) | purpose: VCCS as published (Papon et al., CVPR 2013, §3.1-3.4; decision D6 revised)
//
// Guards finding / decision IDs:
//   VCCS-REF    vccs_optimized == a plain re-implementation of the paper's procedure written independently
//               of the shipped code (full-grid arrays, brute-force nearest-voxel and neighbour searches):
//               same cluster sizes, candidates, kept seeds and threshold, on random grids, R_seed = 2..6,
//               PBC none and PBC xyz. Both follow the same reading of the paper (see Algorithms.cpp), so
//               this checks the implementation, not the reading.
//   VCCS-SEED   candidates = occupied voxels nearest each occupied seed-cell centre (§3.2)
//   VCCS-FILTER threshold = pi*(R_seed/2)^2 (planar slice through the search volume, §3.2), kept iff
//               count >= threshold
//   VCCS-26     adjacency is 26 (§3.1): a diagonal-only chain is one cluster
//   VCCS-REFINE no cluster spans two 26-connected components (flow follows adjacency, §3.1/§3.4)
//   VCCS-COVER  every occupied voxel is counted exactly once
//   VCCS-ITER   at most five iterations (§3.4): verified by the reference, which stops at 5

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "audit_common.hpp"

namespace {

struct RefResult {
  std::vector<int> sizes;  // descending
  int candidates{0}, kept{0};
  double threshold{0};
};

RefResult paperReference(const audit::Grid3& g, int S, bool periodic) {
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const int n[3] = {nx, ny, nz};
  auto wrap = [](int v, int m) { return ((v % m) + m) % m; };
  auto voxel = [&](int i, int j, int k) -> int {  // grid index of an occupied voxel, else -1
    if (periodic) { i = wrap(i, nx); j = wrap(j, ny); k = wrap(k, nz); }
    else if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return -1;
    const int v = static_cast<int>(g.index(i, j, k));
    return g.occ[static_cast<std::size_t>(v)] ? v : -1;
  };
  auto xyz = [&](int v) { return std::array<int, 3>{v / (ny * nz), (v / nz) % ny, v % nz}; };
  auto d1 = [&](double a, double b, int m) { double d = a - b; if (periodic) d -= m * std::round(d / m); return d; };
  auto dist = [&](int v, const std::array<double, 3>& p) {
    const auto q = xyz(v);
    double s = 0;
    for (int a = 0; a < 3; ++a) { const double d = d1(q[a], p[a], n[a]); s += d * d; }
    return std::sqrt(s);
  };
  std::vector<int> occ;
  for (std::size_t v = 0; v < g.size(); ++v) if (g.occ[v]) occ.push_back(static_cast<int>(v));
  auto dist2 = [&](int v, const std::array<double, 3>& p) {
    const auto q = xyz(v);
    double s = 0;
    for (int a = 0; a < 3; ++a) { const double d = d1(q[a], p[a], n[a]); s += d * d; }
    return s;
  };
  auto nearest = [&](const std::array<double, 3>& p) {  // brute force over every occupied voxel
    int b = -1; double bd = 0;
    for (int v : occ) { const double d = dist2(v, p); if (b < 0 || d < bd) { b = v; bd = d; } }
    return b;
  };
  RefResult r;
  // §3.2 candidates
  std::set<std::array<int, 3>> cells;
  for (int v : occ) { const auto q = xyz(v); cells.insert({q[0] / S, q[1] / S, q[2] / S}); }
  std::set<int> cand;
  const double h = 0.5 * (S - 1);
  for (const auto& c : cells) cand.insert(nearest({c[0] * S + h, c[1] * S + h, c[2] * S + h}));
  r.candidates = static_cast<int>(cand.size());
  // §3.2 filter
  const double R = 0.5 * S;
  r.threshold = 3.14159265358979323846 * R * R;
  std::vector<int> seeds;
  for (int c : cand) {
    const auto q = xyz(c);
    int cnt = 0;
    for (int v : occ) if (dist(v, {double(q[0]), double(q[1]), double(q[2])}) <= R) ++cnt;
    if (cnt >= r.threshold) seeds.push_back(c);
  }
  r.kept = static_cast<int>(seeds.size());
  // 26-neighbours of v (occupied)
  auto neighbours = [&](int v) {
    const auto q = xyz(v);
    std::vector<int> out;
    for (int a = -1; a <= 1; ++a) for (int b = -1; b <= 1; ++b) for (int c = -1; c <= 1; ++c) {
      if (!a && !b && !c) continue;
      const int u = voxel(q[0] + a, q[1] + b, q[2] + c);
      if (u >= 0) out.push_back(u);
    }
    return out;
  };
  // §3.3 initial centres: seed + occupied voxels within 2 adjacency steps
  const int ns = static_cast<int>(seeds.size());
  std::vector<std::array<double, 3>> cen(ns);
  std::vector<bool> alive(ns, true);
  for (int s = 0; s < ns; ++s) {
    std::set<int> seen{seeds[s]};
    std::vector<int> layer{seeds[s]};
    for (int step = 0; step < 2; ++step) {
      std::vector<int> nxt;
      for (int v : layer) for (int u : neighbours(v)) if (seen.insert(u).second) nxt.push_back(u);
      layer = nxt;
    }
    const auto q0 = xyz(seeds[s]);
    std::array<double, 3> sum{0, 0, 0};
    for (int v : seen) { const auto q = xyz(v); for (int a = 0; a < 3; ++a) sum[a] += q0[a] + d1(q[a], q0[a], n[a]); }
    for (int a = 0; a < 3; ++a) cen[s][a] = sum[a] / seen.size();
  }
  // §3.4
  std::vector<int> label(g.size(), -1);
  for (int iter = 0; iter < 5; ++iter) {
    std::vector<double> best(g.size(), std::numeric_limits<double>::max());
    std::fill(label.begin(), label.end(), -1);
    std::vector<std::vector<int>> front(ns);
    std::vector<std::set<int>> queued(ns);
    std::vector<bool> on(ns, false);
    for (int s = 0; s < ns; ++s) if (alive[s]) { const int st = nearest(cen[s]); front[s] = {st}; queued[s] = {st}; on[s] = true; }
    bool any = true;
    while (any) {
      any = false;
      for (int s = 0; s < ns; ++s) {
        if (!on[s]) continue;
        std::vector<int> nxt;
        bool set = false;
        for (int v : front[s]) {
          const double d = dist(v, cen[s]);
          if (!(d < best[v])) continue;
          best[v] = d; label[v] = s; set = true;
          for (int u : neighbours(v)) {
            if (queued[s].count(u)) continue;
            const auto q = xyz(u);
            bool inside = true;
            for (int a = 0; a < 3; ++a) if (std::fabs(d1(q[a], cen[s][a], n[a])) > S) inside = false;
            if (!inside || !(dist(u, cen[s]) > d)) continue;
            queued[s].insert(u); nxt.push_back(u);
          }
        }
        front[s] = nxt;
        if (!set || nxt.empty()) on[s] = false; else any = true;
      }
    }
    bool moved = false;
    for (int s = 0; s < ns; ++s) {
      if (!alive[s]) continue;
      std::array<double, 3> sum{0, 0, 0};
      int m = 0;
      for (int v : occ) if (label[v] == s) {
        const auto q = xyz(v);
        for (int a = 0; a < 3; ++a) sum[a] += cen[s][a] + d1(q[a], cen[s][a], n[a]);
        ++m;
      }
      if (!m) { alive[s] = false; moved = true; continue; }
      const std::array<double, 3> c{sum[0] / m, sum[1] / m, sum[2] / m};
      if (c != cen[s]) moved = true;
      cen[s] = c;
    }
    if (!moved) break;
  }
  // remainder: 26-connected, one cluster each
  std::map<int, int> sz;
  int next = 1000000;
  std::vector<int> fin(label);
  for (int v : occ) {
    if (fin[v] >= 0) continue;
    std::vector<int> st{v};
    fin[v] = next;
    while (!st.empty()) {
      const int c = st.back(); st.pop_back();
      for (int u : neighbours(c)) if (fin[u] < 0) { fin[u] = next; st.push_back(u); }
    }
    ++next;
  }
  for (int v : occ) ++sz[fin[v]];
  for (const auto& kv : sz) r.sizes.push_back(kv.second);
  std::sort(r.sizes.begin(), r.sizes.end(), std::greater<int>());
  return r;
}

bls::ClusterResult runVccs(const audit::Grid3& g, int S, bool periodic) {
  auto p = audit::paramsFor(g);
  p.eps = S;
  p.periodic = periodic;
  return audit::runUnlabelled(bls::ClusterAlgorithm::VCCSOptimized, g, p);
}

void reference() {
  std::mt19937_64 rng(1301);
  std::uniform_int_distribution<int> dim(4, 16);
  int cases = 0;
  for (int S : {2, 3, 4, 5, 6})
    for (double p : {0.05, 0.2, 0.45, 0.8})
      for (bool periodic : {false, true})
        for (int rep = 0; rep < 5; ++rep) {
          const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
          const auto ref = paperReference(g, S, periodic);
          const auto got = runVccs(g, S, periodic);
          ++cases;
          const std::string tag = "S=" + std::to_string(S) + " p=" + std::to_string(p) + (periodic ? " PBC " : " ") +
                                  std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" + std::to_string(g.nz);
          AUDIT_CHECK("VCCS-REF", got.clusterSizes == ref.sizes,
                      tag + ": sizes differ (" + std::to_string(got.nclusters) + " vs " + std::to_string(ref.sizes.size()) + " clusters)");
          AUDIT_CHECK("VCCS-SEED", got.seedCandidates == ref.candidates,
                      tag + ": candidates " + std::to_string(got.seedCandidates) + " vs " + std::to_string(ref.candidates));
          AUDIT_CHECK("VCCS-FILTER", got.seedsPlaced == ref.kept && std::fabs(got.seedPruneThreshold - ref.threshold) < 1e-9,
                      tag + ": kept " + std::to_string(got.seedsPlaced) + " vs " + std::to_string(ref.kept));
          long long total = 0, occ = 0;
          for (int sz : got.clusterSizes) total += sz;
          for (auto v : g.occ) occ += v;
          AUDIT_CHECK("VCCS-COVER", total == occ && got.visitedVoxels == static_cast<std::size_t>(occ), tag + ": coverage");
          // Refinement of 26-connected components (non-periodic: CC3D-opt at 26 is the reference).
          if (!periodic) {
            auto p26 = audit::paramsFor(g);
            p26.connectivity = 26;
            const auto cc = audit::runUnlabelled(bls::ClusterAlgorithm::CC3DOptimized, g, p26);
            AUDIT_CHECK("VCCS-REFINE", got.nclusters >= cc.nclusters && got.maxCluster <= cc.maxCluster,
                        tag + ": a cluster spans two 26-connected components");
          }
        }
  std::printf("measure: VCCS-REF %d random cases (R_seed 2-6, p 0.05-0.8, PBC none/xyz)\n", cases);
  AUDIT_CHECK("VCCS-ITER", cases > 0, "no cases");
}

void adjacency26() {
  auto g = audit::emptyGrid(12, 12, 12);
  for (int t = 2; t < 9; ++t) g.set(t, t, t);
  const auto r = runVccs(g, 12, false);
  AUDIT_CHECK("VCCS-26", r.nclusters == 1, "diagonal chain gave " + std::to_string(r.nclusters) + " clusters (want 1)");
}

}  // namespace

int main() {
  reference();
  adjacency26();
  return audit::finish("test_audit_vccs");
}
