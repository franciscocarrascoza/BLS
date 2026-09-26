// Created: 2026-09-26T17:54+02:00 | by: CC audit (manager) | purpose: WP3 tests — periodic adjacency (deck PBC xyz) in the exact comparison methods
//
// Guards decision / claim IDs:
//   PBC-M      traditional_dfs, gcbd, cc3d_optimized (6 and 26) and rle_ccl_optimized under PBC xyz label
//              exactly the periodic components (union-find oracle with wrap-around), on random grids that
//              include 1-, 2- and 3-voxel dimensions; skip_dfs at unit stride (Tier 1, D8) gives the same
//              component count and sizes
//   PBC-TRANS  periodic translation covariance: rolling the grid by any vector leaves every method's
//              component sizes unchanged (independent of the oracle)
//   PBC-NONE   periodic = false still gives the non-periodic partition (the PBC none instantiation)
//   PBC-REFUSE every method without periodic adjacency throws under PBC instead of returning a
//              non-periodic answer: withdrawn cc3d/rle_ccl/vccs, dbscan/hdbscan (D4/D5), k-means (D11),
//              hierarchical, vccs_optimized (until WP4)

#include <algorithm>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "audit_common.hpp"

namespace {

using bls::ClusterAlgorithm;

audit::Grid3 roll(const audit::Grid3& g, int ti, int tj, int tk) {
  audit::Grid3 r = audit::emptyGrid(g.nx, g.ny, g.nz);
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k)
        if (g.at(i, j, k)) r.set((i + ti) % g.nx, (j + tj) % g.ny, (k + tk) % g.nz);
  return r;
}

std::vector<int> sortedSizes(const std::vector<int>& labels) { return audit::sizesOf(labels); }

struct Labelled {
  ClusterAlgorithm a;
  int connectivity;
  const char* name;
};
const Labelled kLabelled[] = {
    {ClusterAlgorithm::TraditionalDFS, 6, "traditional_dfs"},
    {ClusterAlgorithm::GCBD, 6, "gcbd"},
    {ClusterAlgorithm::CC3DOptimized, 6, "cc3d_optimized/6"},
    {ClusterAlgorithm::CC3DOptimized, 26, "cc3d_optimized/26"},
    {ClusterAlgorithm::RLECCLOptimized, 6, "rle_ccl_optimized"},
};

void periodicPartitions() {
  std::mt19937_64 rng(20260926);
  std::uniform_int_distribution<int> dim(1, 16), shift(0, 15);
  int grids = 0;
  for (double p : {0.1, 0.3, 0.5, 0.7}) {
    for (int rep = 0; rep < 60; ++rep) {
      const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
      ++grids;
      const auto oracle6 = audit::canonical(audit::periodicComponents(g, 6));
      const auto oracle26 = audit::canonical(audit::periodicComponents(g, 26));
      const auto r = roll(g, shift(rng) % g.nx, shift(rng) % g.ny, shift(rng) % g.nz);
      const std::string shape = std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" + std::to_string(g.nz) +
                                " p=" + std::to_string(p);
      for (const auto& m : kLabelled) {
        auto params = audit::paramsFor(g);
        params.connectivity = m.connectivity;
        params.periodic = true;
        std::vector<int> labels, rolled, open;
        const auto res = audit::runLabelled(m.a, g, labels, params);
        AUDIT_CHECK("PBC-M", audit::labelsCoverOccupancy(g, labels) &&
                                 audit::canonical(labels) == (m.connectivity == 26 ? oracle26 : oracle6) &&
                                 res.nclusters == static_cast<int>(audit::sizesOf(labels).size()),
                    std::string(m.name) + " " + shape + ": periodic partition != oracle (nclusters " +
                        std::to_string(res.nclusters) + ")");
        audit::runLabelled(m.a, r, rolled, params);
        AUDIT_CHECK("PBC-TRANS", sortedSizes(rolled) == sortedSizes(labels),
                    std::string(m.name) + " " + shape + ": sizes change under a periodic roll");
        params.periodic = false;
        audit::runLabelled(m.a, g, open, params);
        std::vector<int> dfsOpen;
        auto p6 = audit::paramsFor(g);
        p6.connectivity = m.connectivity;
        audit::runLabelled(m.connectivity == 26 ? ClusterAlgorithm::CC3DOptimized : ClusterAlgorithm::TraditionalDFS, g,
                           dfsOpen, p6);
        AUDIT_CHECK("PBC-NONE", audit::canonical(open) == audit::canonical(dfsOpen),
                    std::string(m.name) + " " + shape + ": PBC none partition changed");
      }
      // skip_dfs at unit stride: no label output; count and sizes.
      auto params = audit::paramsFor(g);
      params.skipDfsJumpDistance = 1;
      params.periodic = true;
      const auto s = audit::runUnlabelled(ClusterAlgorithm::SkipDFS, g, params);
      const auto want = audit::sizesOf(audit::periodicComponents(g, 6));
      AUDIT_CHECK("PBC-M", s.nclusters == static_cast<int>(want.size()) && s.clusterSizes == want,
                  "skip_dfs s=1 " + shape + ": " + std::to_string(s.nclusters) + " clusters vs oracle " +
                      std::to_string(want.size()));
      const auto sr = audit::runUnlabelled(ClusterAlgorithm::SkipDFS, r, params);
      AUDIT_CHECK("PBC-TRANS", sr.clusterSizes == s.clusterSizes, "skip_dfs s=1 " + shape + ": sizes change under a roll");
    }
  }
  std::printf("measure: PBC-M %d random periodic grids (dims 1-16, p 0.1-0.7) x 5 labelled methods + skip_dfs\n", grids);
}

void refusals() {
  std::mt19937_64 rng(7);
  const auto g = audit::randomGrid(rng, 8, 8, 8, 0.3);
  const ClusterAlgorithm refused[] = {ClusterAlgorithm::CC3D,     ClusterAlgorithm::RLECCL,
                                      ClusterAlgorithm::VCCS,     ClusterAlgorithm::VCCSOptimized,
                                      ClusterAlgorithm::DBSCAN,   ClusterAlgorithm::HDBSCAN,
                                      ClusterAlgorithm::KMeans,   ClusterAlgorithm::Hierarchical};
  for (auto a : refused) {
    auto params = audit::paramsFor(g);
    params.periodic = true;
    params.k = 3;
    bool threw = false;
    try {
      audit::runUnlabelled(a, g, params);
    } catch (const std::runtime_error& e) {
      threw = std::string(e.what()).find("PBC is not implemented") != std::string::npos;
    }
    AUDIT_CHECK("PBC-REFUSE", threw && !bls::supportsPeriodic(a),
                bls::algorithmToString(a) + " did not refuse PBC xyz");
  }
  for (auto a : {ClusterAlgorithm::BLS, ClusterAlgorithm::TraditionalDFS, ClusterAlgorithm::SkipDFS,
                 ClusterAlgorithm::GCBD, ClusterAlgorithm::CC3DOptimized, ClusterAlgorithm::RLECCLOptimized}) {
    AUDIT_CHECK("PBC-REFUSE", bls::supportsPeriodic(a), bls::algorithmToString(a) + " should support PBC xyz");
  }
}

}  // namespace

int main() {
  periodicPartitions();
  refusals();
  return audit::finish("test_audit_pbc_methods");
}
