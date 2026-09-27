// Created: 2026-09-27T18:06+02:00 | by: CC audit (manager) | purpose: HDBSCAN* as published (Campello, Moulavi, Sander, PAKDD 2013) — decision D4/D5 revised
//
// Guards finding / decision IDs:
//   HDB-REF   hdbscan == a plain reference that builds each hierarchy level directly as the DBSCAN*
//             partition on the COMPLETE mutual-reachability graph (the paper's Proposition 1), not from
//             an MST: objects with d_core < w, joined by d_mreach < w, for every distinct level w (all
//             ties at once), then Algorithm 2 (m_clSize), Eq. (3) and Algorithm 3 (root excluded).
//             Random grids, m_pts 2..5, m_clSize 1..5, PBC none/xyz: same cluster count and sizes.
//   HDB-BLOBS well-separated blobs are separated (the old code returned one cluster: finding S1-9)
//   HDB-ROOT  a single blob that never splits yields no cluster (root excluded, Algorithm 3)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "audit_common.hpp"

namespace {

std::vector<int> referenceHdbscan(const audit::Grid3& g, int mpts, int mcl, bool periodic) {
  std::vector<std::array<int, 3>> P;
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k)
        if (g.at(i, j, k)) P.push_back({i, j, k});
  const int n = static_cast<int>(P.size());
  if (n < mpts) return {};
  const int dim[3] = {g.nx, g.ny, g.nz};
  auto d2 = [&](int a, int b) {
    long long s = 0;
    for (int t = 0; t < 3; ++t) {
      int d = std::abs(P[a][t] - P[b][t]);
      if (periodic) d = std::min(d, dim[t] - d);
      s += 1LL * d * d;
    }
    return s;
  };
  std::vector<long long> core(n);
  for (int a = 0; a < n; ++a) {
    std::vector<long long> all;
    for (int b = 0; b < n; ++b) all.push_back(d2(a, b));
    std::sort(all.begin(), all.end());
    core[a] = all[mpts - 1];
  }
  auto mr = [&](int a, int b) { return std::max({core[a], core[b], d2(a, b)}); };
  std::vector<long long> levels(core.begin(), core.end());
  for (int a = 0; a < n; ++a)
    for (int b = a + 1; b < n; ++b) levels.push_back(mr(a, b));
  std::sort(levels.begin(), levels.end(), std::greater<long long>());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

  struct Cl { int parent; double lb; double S; std::vector<int> kids; int birthSize; };
  std::vector<Cl> cl{{-1, 0.0, 0.0, {}, n}};
  std::vector<int> of(n, 0);
  for (long long w : levels) {
    const double lam = 1.0 / std::sqrt(static_cast<double>(w));
    // DBSCAN* partition after removing every edge of weight >= w (Proposition 1, complete graph).
    std::vector<int> comp(n, -1);
    int nc = 0;
    for (int a = 0; a < n; ++a) {
      if (comp[a] >= 0) continue;
      comp[a] = nc;
      std::vector<int> st{a};
      while (!st.empty()) {
        int x = st.back(); st.pop_back();
        for (int b = 0; b < n; ++b)
          if (comp[b] < 0 && mr(x, b) < w) { comp[b] = nc; st.push_back(b); }
      }
      ++nc;
    }
    std::map<int, std::vector<int>> byCluster;  // current cluster -> its members
    for (int a = 0; a < n; ++a) if (of[a] >= 0) byCluster[of[a]].push_back(a);
    for (auto& [c, mem] : byCluster) {
      std::map<int, std::vector<int>> subs;
      for (int a : mem) subs[comp[a]].push_back(a);
      std::vector<std::vector<int>> keep;
      for (auto& [k, sub] : subs) {
        const bool hasEdge = sub.size() >= 2 || core[sub[0]] < w;
        if (!(static_cast<int>(sub.size()) < mcl || (mcl == 1 && !hasEdge))) keep.push_back(sub);
      }
      if (keep.size() == 1 && keep[0].size() == mem.size()) continue;  // unchanged at this level
      const double gain = lam - cl[c].lb;
      if (keep.size() == 1) {
        cl[c].S += static_cast<double>(mem.size() - keep[0].size()) * gain;
        for (int a : mem) of[a] = -1;
        for (int a : keep[0]) of[a] = c;
      } else {
        cl[c].S += static_cast<double>(mem.size()) * gain;
        for (int a : mem) of[a] = -1;
        for (auto& sub : keep) {
          const int id = static_cast<int>(cl.size());
          cl.push_back({c, lam, 0.0, {}, static_cast<int>(sub.size())});
          cl[c].kids.push_back(id);
          for (int a : sub) of[a] = id;
        }
      }
    }
  }
  const int nc = static_cast<int>(cl.size());
  std::vector<double> sh(nc, 0.0);
  std::vector<bool> keepC(nc, false);
  for (int c = nc - 1; c >= 1; --c) {
    double sum = 0;
    for (int k : cl[c].kids) sum += sh[k];
    if (cl[c].kids.empty() || !(cl[c].S < sum)) { sh[c] = cl[c].S; keepC[c] = true; }
    else { sh[c] = sum; keepC[c] = false; }
  }
  std::vector<int> sizes;
  std::function<void(int)> pick = [&](int c) {
    if (keepC[c]) { sizes.push_back(cl[c].birthSize); return; }
    for (int k : cl[c].kids) pick(k);
  };
  for (int k : cl[0].kids) pick(k);
  std::sort(sizes.begin(), sizes.end(), std::greater<int>());
  return sizes;
}

bls::ClusterResult run(const audit::Grid3& g, int mpts, int mcl, bool periodic) {
  auto p = audit::paramsFor(g);
  p.minSamples = mpts;
  p.minClusterSize = mcl;
  p.periodic = periodic;
  return audit::runUnlabelled(bls::ClusterAlgorithm::HDBSCAN, g, p);
}

void reference() {
  std::mt19937_64 rng(2013);
  std::uniform_int_distribution<int> dim(3, 9);
  int cases = 0, withClusters = 0;
  for (int mpts : {2, 3, 5})
    for (int mcl : {1, 2, 3, 5})
      for (double p : {0.15, 0.3, 0.5})
        for (bool periodic : {false, true})
          for (int rep = 0; rep < 3; ++rep) {
            const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
            const auto want = referenceHdbscan(g, mpts, mcl, periodic);
            const auto got = run(g, mpts, mcl, periodic);
            ++cases;
            if (!want.empty()) ++withClusters;
            AUDIT_CHECK("HDB-REF", got.clusterSizes == want,
                        "m_pts=" + std::to_string(mpts) + " m_clSize=" + std::to_string(mcl) + " p=" + std::to_string(p) +
                            (periodic ? " PBC " : " ") + std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" +
                            std::to_string(g.nz) + ": " + std::to_string(got.nclusters) + " clusters vs reference " +
                            std::to_string(want.size()));
          }
  std::printf("measure: HDB-REF %d random cases (%d with at least one cluster)\n", cases, withClusters);
}

void blobs() {
  auto g = audit::emptyGrid(64, 8, 8);
  for (int b = 0; b < 5; ++b)
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) g.set(2 + 12 * b + i, 2 + j, 2 + k);
  const auto r = run(g, 5, 5, false);
  AUDIT_CHECK("HDB-BLOBS", r.nclusters == 5 && r.maxCluster == 27, "5 blobs gave " + std::to_string(r.nclusters));
  auto one = audit::emptyGrid(8, 8, 8);
  for (int i = 2; i < 5; ++i) for (int j = 2; j < 5; ++j) for (int k = 2; k < 5; ++k) one.set(i, j, k);
  AUDIT_CHECK("HDB-ROOT", run(one, 5, 5, false).nclusters == 0, "a single blob should give no cluster (root excluded)");
}

}  // namespace

int main() {
  reference();
  blobs();
  return audit::finish("test_audit_hdbscan");
}
