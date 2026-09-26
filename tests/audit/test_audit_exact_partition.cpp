// Created: 2026-09-26T12:47+02:00 | by: CC audit (manager) | purpose: exact-tier partition equality on 10^4 random + structured grids
//
// Guards claims-ledger IDs:
//   F3   "any method claiming exactness must reproduce their output component for component" (§1.1 l.117-118)
//   C-2  exact tier returns identical counts and maximum sizes (§4.2 l.680-683, Table 4)
//   S-2  unit tests cover cross-algorithm output agreement (§3.1 l.393-397)
//
// Every exact method's labelling is compared with traditional_dfs by CANONICAL
// PARTITION (min-voxel relabelling), not by count. Also checked per method: the
// returned nclusters / maxCluster / clusterSizes agree with its own labels, and
// labels cover exactly the occupied voxels.
//
// Reported methods: cc3d_optimized (CC3D), rle_ccl_optimized (RLE-CCL), gcbd.
// Withdrawn methods cc3d, rle_ccl are checked under a separate tag (WITHDRAWN) so a
// failure there cannot be confused with a reported-method failure.
//
// Usage: test_audit_exact_partition [n_random_grids]   (default 10000)

#include <cstdlib>
#include <string>

#include "audit_common.hpp"

namespace {

using audit::Grid3;

struct Method {
  bls::ClusterAlgorithm algo;
  const char* name;
  const char* claim;  // ledger tag this method's checks count toward
};

const Method kMethods[] = {
    {bls::ClusterAlgorithm::CC3DOptimized, "cc3d_optimized", "F3/C-2"},
    {bls::ClusterAlgorithm::RLECCLOptimized, "rle_ccl_optimized", "F3/C-2"},
    {bls::ClusterAlgorithm::GCBD, "gcbd", "F3/C-2"},
    {bls::ClusterAlgorithm::CC3D, "cc3d(withdrawn)", "WITHDRAWN"},
    {bls::ClusterAlgorithm::RLECCL, "rle_ccl(withdrawn)", "WITHDRAWN"},
};

std::string dims(const Grid3& g) {
  return std::to_string(g.nx) + "x" + std::to_string(g.ny) + "x" + std::to_string(g.nz);
}

void compareAll(const Grid3& g, const std::string& caseName) {
  std::vector<int> ref;
  bls::ClusterResult rr = audit::runLabelled(bls::ClusterAlgorithm::TraditionalDFS, g, ref);
  AUDIT_CHECK("S-2", audit::labelsCoverOccupancy(g, ref), "dfs labels do not cover occupancy on " + caseName);
  const auto refCanon = audit::canonical(ref);
  const auto refSizes = audit::sizesOf(ref);
  AUDIT_CHECK("S-2", rr.nclusters == static_cast<int>(refSizes.size()),
              "dfs nclusters != distinct labels on " + caseName);

  for (const Method& m : kMethods) {
    std::vector<int> lab;
    bls::ClusterResult r = audit::runLabelled(m.algo, g, lab);
    const std::string where = std::string(m.name) + " on " + caseName + " (" + dims(g) + ")";
    AUDIT_CHECK(m.claim, audit::labelsCoverOccupancy(g, lab), "labels do not cover occupancy: " + where);
    AUDIT_CHECK(m.claim, audit::canonical(lab) == refCanon, "partition differs from DFS: " + where);
    std::vector<int> sz = r.clusterSizes;
    std::sort(sz.begin(), sz.end(), std::greater<int>());
    AUDIT_CHECK(m.claim, r.nclusters == rr.nclusters, "nclusters differs from DFS: " + where);
    AUDIT_CHECK(m.claim, sz == refSizes, "size list differs from DFS: " + where);
    AUDIT_CHECK(m.claim, r.maxCluster == (refSizes.empty() ? 0 : refSizes.front()),
                "maxCluster differs from DFS: " + where);
    AUDIT_CHECK(m.claim, audit::sizesOf(lab) == sz, "reported sizes disagree with own labels: " + where);
  }
}

// Structured cases: degenerate shapes and components touching every face.
void structured() {
  // 1x1x1 empty and full
  for (int v = 0; v <= 1; ++v) {
    Grid3 g = audit::emptyGrid(1, 1, 1);
    g.occ[0] = static_cast<uint8_t>(v);
    compareAll(g, "1x1x1 v=" + std::to_string(v));
  }
  // 1xNxN, Nx1xN, NxNx1 random slabs
  std::mt19937_64 rng(20260926);
  for (int n : {2, 7, 31}) {
    for (double p : {0.3, 0.6}) {
      compareAll(audit::randomGrid(rng, 1, n, n, p), "1xNxN");
      compareAll(audit::randomGrid(rng, n, 1, n, p), "Nx1xN");
      compareAll(audit::randomGrid(rng, n, n, 1, p), "NxNx1");
      compareAll(audit::randomGrid(rng, 1, 1, n, p), "1x1xN");
    }
  }
  // full grid; hollow shell touching all six faces; 3D checkerboard (only diagonal contacts)
  for (int n : {1, 2, 5, 16}) {
    Grid3 full = audit::emptyGrid(n, n + 1, n + 2);
    std::fill(full.occ.begin(), full.occ.end(), 1);
    compareAll(full, "full");
    Grid3 shell = audit::emptyGrid(n + 2, n + 3, n + 4);
    for (int i = 0; i < shell.nx; ++i)
      for (int j = 0; j < shell.ny; ++j)
        for (int k = 0; k < shell.nz; ++k)
          if (i == 0 || j == 0 || k == 0 || i == shell.nx - 1 || j == shell.ny - 1 || k == shell.nz - 1)
            shell.set(i, j, k);
    compareAll(shell, "shell");
    Grid3 cb = audit::emptyGrid(n + 3, n + 3, n + 3);
    for (int i = 0; i < cb.nx; ++i)
      for (int j = 0; j < cb.ny; ++j)
        for (int k = 0; k < cb.nz; ++k)
          if ((i + j + k) % 2 == 0) cb.set(i, j, k);
    compareAll(cb, "checkerboard");
  }
  // U-shapes that merge late in raster order (classic two-pass stress), in all 3 orientations
  for (int axis = 0; axis < 3; ++axis) {
    Grid3 g = audit::emptyGrid(9, 9, 9);
    for (int a = 0; a < 9; ++a)
      for (int b = 0; b < 9; ++b) {
        const bool on = (b == 0) || (a == 0) || (a == 8) || (b == 4 && a % 2 == 0);
        if (!on) continue;
        if (axis == 0) g.set(4, a, b);
        if (axis == 1) g.set(a, 4, b);
        if (axis == 2) g.set(a, b, 4);
      }
    compareAll(g, "U-shape axis " + std::to_string(axis));
  }
}

}  // namespace

int main(int argc, char** argv) {
  const int nRandom = argc > 1 ? std::atoi(argv[1]) : 10000;
  structured();
  std::mt19937_64 rng(0xB15A0D17ULL);  // fixed seed
  std::uniform_int_distribution<int> dim(1, 24);
  std::uniform_real_distribution<double> occ(0.05, 0.6);
  for (int t = 0; t < nRandom; ++t) {
    const int nx = dim(rng), ny = dim(rng), nz = dim(rng);
    const double p = occ(rng);
    compareAll(audit::randomGrid(rng, nx, ny, nz, p), "random#" + std::to_string(t));
  }
  std::printf("random grids: %d (seed 0xB15A0D17, dims 1..24, occupancy 0.05..0.6)\n", nRandom);
  return audit::finish("test_audit_exact_partition");
}
