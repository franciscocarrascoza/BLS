// Created: 2026-09-26T12:48+02:00 | by: CC audit (manager) | purpose: reduction tests (degenerate parameters must reproduce N6 CCL) + published-form sanity checks
//
// Guards claims-ledger IDs:
//   M-9/C-6  Skip-DFS at s = 1 is N6 DFS (§2.3 l.311-316; §4.2 l.691 "unit stride")        -> partition (count+sizes)
//   F5       DBSCAN, eps in [1, sqrt2) voxel, minPts = 1  == N6 CCL                         -> count + size multiset
//   F6       single-linkage hierarchical, threshold in [1, sqrt2) == N6 CCL                 -> count + size multiset
//   M-8      BLS with dNN well below one voxel (every voxel probed) == DFS                  -> canonical partition
//   F5-PUB   DBSCAN follows Ester et al. 1996 (expansion from core points only)             -> vs textbook reference
//   C-8-PUB  HDBSCAN separates well-separated blobs (R3: single-component output)          -> 5 blobs must give 5
//
// DBSCAN, hierarchical and Skip-DFS provide no per-voxel labels in this code base
// (supportsLabels() == false), so for them the check is count + full sorted size
// multiset. That is weaker than partition equality; it is stated as such.

#include <cmath>
#include <cstdlib>
#include <string>

#include "audit_common.hpp"

namespace {

using audit::Grid3;

bls::ClusterResult dfsOf(const Grid3& g) { return audit::runUnlabelled(bls::ClusterAlgorithm::TraditionalDFS, g, audit::paramsFor(g)); }

bool sameCountAndSizes(const bls::ClusterResult& a, const bls::ClusterResult& b) {
  return a.nclusters == b.nclusters && a.clusterSizes == b.clusterSizes;
}

std::string desc(const Grid3& g, const std::string& what) {
  return what + " on " + std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" + std::to_string(g.nz);
}

// ---------------------------------------------------------------------------
// Textbook DBSCAN (Ester, Kriegel, Sander, Xu 1996), O(N^2), voxel-index Euclidean
// distance, |N_eps(p)| counts p itself. Points are visited in the same x-major
// scan order as the library, so border points reachable from two clusters are
// assigned identically. Returns descending cluster sizes (noise excluded).
// ---------------------------------------------------------------------------
std::vector<int> referenceDbscan(const Grid3& g, double eps, int minPts) {
  struct P { int i, j, k; };
  std::vector<P> pts;
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k)
        if (g.at(i, j, k)) pts.push_back({i, j, k});
  const int n = static_cast<int>(pts.size());
  auto region = [&](int p) {
    std::vector<int> out;
    for (int q = 0; q < n; ++q) {
      const double d = std::sqrt(double((pts[p].i - pts[q].i) * (pts[p].i - pts[q].i) +
                                        (pts[p].j - pts[q].j) * (pts[p].j - pts[q].j) +
                                        (pts[p].k - pts[q].k) * (pts[p].k - pts[q].k)));
      if (d <= eps) out.push_back(q);
    }
    return out;
  };
  const int UNDEF = -2, NOISE = -1;
  std::vector<int> lab(n, UNDEF);
  int c = 0;
  for (int p = 0; p < n; ++p) {
    if (lab[p] != UNDEF) continue;
    std::vector<int> nb = region(p);
    if (static_cast<int>(nb.size()) < minPts) { lab[p] = NOISE; continue; }
    lab[p] = c;
    std::vector<int> seeds(nb.begin(), nb.end());
    for (std::size_t s = 0; s < seeds.size(); ++s) {
      const int q = seeds[s];
      if (lab[q] == NOISE) lab[q] = c;  // border point
      if (lab[q] != UNDEF) continue;
      lab[q] = c;
      std::vector<int> nq = region(q);
      if (static_cast<int>(nq.size()) >= minPts) seeds.insert(seeds.end(), nq.begin(), nq.end());
    }
    ++c;
  }
  std::vector<int> sizes(c, 0);
  for (int l : lab) if (l >= 0) ++sizes[l];
  std::sort(sizes.begin(), sizes.end(), std::greater<int>());
  return sizes;
}

void skipDfsUnitStride(std::mt19937_64& rng, int n) {
  int differAtS3 = 0;
  for (int t = 0; t < n; ++t) {
    Grid3 g = audit::randomGrid(rng, 4 + t % 20, 4 + (t * 7) % 20, 4 + (t * 13) % 20, 0.05 + 0.5 * (t % 11) / 10.0);
    const auto ref = dfsOf(g);
    bls::ClusterParams p = audit::paramsFor(g);
    p.skipDfsJumpDistance = 1;
    AUDIT_CHECK("M-9/C-6", sameCountAndSizes(audit::runUnlabelled(bls::ClusterAlgorithm::SkipDFS, g, p), ref),
                desc(g, "skip_dfs s=1 != DFS"));
    p.skipDfsJumpDistance = 3;
    if (!sameCountAndSizes(audit::runUnlabelled(bls::ClusterAlgorithm::SkipDFS, g, p), ref)) ++differAtS3;
  }
  std::printf("info: skip_dfs s=3 (the stride every campaign used) differs from DFS on %d of %d grids\n", differAtS3, n);
}

void dbscanReduction(std::mt19937_64& rng, int n) {
  for (double eps : {1.0, 1.2, 1.41}) {
    for (int t = 0; t < n; ++t) {
      Grid3 g = audit::randomGrid(rng, 3 + t % 18, 3 + (t * 5) % 18, 3 + (t * 11) % 18, 0.05 + 0.55 * (t % 7) / 6.0);
      bls::ClusterParams p = audit::paramsFor(g);
      p.eps = eps;
      p.minPts = 1;
      AUDIT_CHECK("F5", sameCountAndSizes(audit::runUnlabelled(bls::ClusterAlgorithm::DBSCAN, g, p), dfsOf(g)),
                  desc(g, "dbscan eps=" + std::to_string(eps) + " minPts=1 != N6 CCL"));
    }
  }
}

void hierarchicalReduction(std::mt19937_64& rng, int n) {
  for (double thr : {1.0, 1.3}) {
    for (int t = 0; t < n; ++t) {
      Grid3 g = audit::randomGrid(rng, 2 + t % 12, 2 + (t * 5) % 12, 2 + (t * 7) % 12, 0.05 + 0.5 * (t % 6) / 5.0);
      bls::ClusterParams p = audit::paramsFor(g);
      p.threshold = thr;
      AUDIT_CHECK("F6", sameCountAndSizes(audit::runUnlabelled(bls::ClusterAlgorithm::Hierarchical, g, p), dfsOf(g)),
                  desc(g, "hierarchical threshold=" + std::to_string(thr) + " != N6 CCL"));
    }
  }
}

struct LatticeCase { const char* name; bls::LatticeType t; bls::CenteringType c; };
const LatticeCase kLattices[] = {
    {"cubic-P", bls::LatticeType::Cubic, bls::CenteringType::P},
    {"cubic-I", bls::LatticeType::Cubic, bls::CenteringType::I},
    {"cubic-F", bls::LatticeType::Cubic, bls::CenteringType::F},
    {"hex-P", bls::LatticeType::Hexagonal, bls::CenteringType::P},
    {"hex-I", bls::LatticeType::Hexagonal, bls::CenteringType::I},
    {"hex-F", bls::LatticeType::Hexagonal, bls::CenteringType::F},
    {"tri-P", bls::LatticeType::Triclinic, bls::CenteringType::P},
    {"tri-I", bls::LatticeType::Triclinic, bls::CenteringType::I},
    {"tri-F", bls::LatticeType::Triclinic, bls::CenteringType::F},
};

// Fraction of voxels of an nx*ny*nz grid that carry a probe (full occupancy).
double probeCoverage(const bls::LatticeDescriptor& d, double dnn, int n) {
  Grid3 full = audit::emptyGrid(n, n, n);
  std::fill(full.occ.begin(), full.occ.end(), 1);
  const bls::Mat3 sb = d.basis * (dnn / d.dmin);
  bls::Enumerator en(sb, d.offsets, n, n, n, full.occ);
  return static_cast<double>(en.count()) / static_cast<double>(full.size());
}

void blsDegenerate(std::mt19937_64& rng, int n) {
  for (const auto& lc : kLattices) {
    const bls::LatticeDescriptor d = audit::latticeFor(lc.t, lc.c);
    for (double dnn : {0.5, 0.7, 1.0}) {
      const double cov = probeCoverage(d, dnn, 16);
      std::printf("info: %-8s dNN=%.2f voxel: probe coverage of a full 16^3 grid = %.4f\n", lc.name, dnn, cov);
      if (cov < 1.0) continue;  // the reduction is only claimed when every voxel is probed
      for (int t = 0; t < n; ++t) {
        Grid3 g = audit::randomGrid(rng, 3 + t % 17, 3 + (t * 3) % 17, 3 + (t * 7) % 17, 0.05 + 0.55 * (t % 5) / 4.0);
        std::vector<int> ref;
        audit::runLabelled(bls::ClusterAlgorithm::TraditionalDFS, g, ref);
        audit::BlsRaw b = audit::runBlsRaw(g, d, dnn);
        AUDIT_CHECK("M-8", audit::canonical(b.labels) == audit::canonical(ref),
                    desc(g, std::string("BLS ") + lc.name + " dNN=" + std::to_string(dnn) + " != DFS partition"));
      }
    }
  }
}

// DBSCAN published form: a 5^3 solid cube with a 4-voxel chain off one face centre,
// eps = 1.0, minPts = 7. Only the 27 interior voxels are core (6 neighbours + self);
// Ester et al. attach border points reachable from a core point and do NOT expand
// from border points.
void dbscanPublishedForm() {
  Grid3 g = audit::emptyGrid(12, 7, 7);
  for (int i = 1; i <= 5; ++i)
    for (int j = 1; j <= 5; ++j)
      for (int k = 1; k <= 5; ++k) g.set(i, j, k);
  for (int c = 6; c <= 9; ++c) g.set(c, 3, 3);  // chain off the x = 5 face centre
  bls::ClusterParams p = audit::paramsFor(g);
  p.eps = 1.0;
  p.minPts = 7;
  const auto got = audit::runUnlabelled(bls::ClusterAlgorithm::DBSCAN, g, p);
  const auto want = referenceDbscan(g, 1.0, 7);
  std::string ws, gs;
  for (int s : want) ws += std::to_string(s) + " ";
  for (int s : got.clusterSizes) gs += std::to_string(s) + " ";
  std::printf("info: DBSCAN cube+chain: textbook sizes [%s], code sizes [%s]\n", ws.c_str(), gs.c_str());
  AUDIT_CHECK("F5-PUB", got.clusterSizes == want, "dbscan differs from Ester et al. 1996 on cube+chain");
  // Also on the campaign parameters (eps 3.0, minPts 10) over small random grids.
  std::mt19937_64 rng(77);
  for (int t = 0; t < 60; ++t) {
    Grid3 r = audit::randomGrid(rng, 8 + t % 6, 8, 8 + (t * 3) % 6, 0.03 + 0.04 * (t % 5));
    bls::ClusterParams q = audit::paramsFor(r);
    q.eps = 3.0;
    q.minPts = 10;
    AUDIT_CHECK("F5-PUB", audit::runUnlabelled(bls::ClusterAlgorithm::DBSCAN, r, q).clusterSizes == referenceDbscan(r, 3.0, 10),
                desc(r, "dbscan(eps 3, minPts 10) differs from textbook DBSCAN"));
  }
}

// HDBSCAN: five solid 3^3 cubes, 12 voxels apart along x, campaign parameters
// (min_cluster_size 5, min_samples 5). Any HDBSCAN* with excess-of-mass selection
// (root excluded) returns 5 clusters here.
void hdbscanBlobs() {
  Grid3 g = audit::emptyGrid(5 * 15, 7, 7);
  for (int b = 0; b < 5; ++b)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) g.set(2 + 15 * b + i, 2 + j, 2 + k);
  bls::ClusterParams p = audit::paramsFor(g);
  p.minClusterSize = 5;
  p.minSamples = 5;
  const auto r = audit::runUnlabelled(bls::ClusterAlgorithm::HDBSCAN, g, p);
  std::printf("info: HDBSCAN 5 separated blobs -> %d cluster(s)\n", r.nclusters);
  AUDIT_CHECK("C-8-PUB", r.nclusters == 5, "hdbscan does not separate 5 well-separated blobs (got " + std::to_string(r.nclusters) + ")");
}

}  // namespace

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 1000;
  std::mt19937_64 rng(0x5EEDF00DULL);
  skipDfsUnitStride(rng, 2 * n);
  dbscanReduction(rng, n / 2);
  hierarchicalReduction(rng, n / 4);
  blsDegenerate(rng, n / 10);
  dbscanPublishedForm();
  hdbscanBlobs();
  return audit::finish("test_audit_reductions");
}
