// Created: 2026-09-28T09:17+02:00 | by: CC audit (manager) | purpose: BLS stage changes of 28.09.26 — separable probe loop (BLS-PERF2), REFINEMENT DFS, stage timers
//
// Guards claims-ledger / finding / decision IDs:
//   SEP-ID    the probe evaluation (separable loop for diagonal bases, merged runs) returns the
//             identical seed list, in the identical order, and the identical probes(), as the
//             row-by-row loop that preceded it (Enumerator::rowWise / rowWisePeriodic), for every
//             lattice x centering, dNN below and above one voxel, occupancies, lattice origins
//             (including half-integer ones, where rounding ties) and grid shapes; periodic and not
//   SEP-ALL   the lattice-only constructor (no occupancy) returns the same sites as the
//             occupancy constructor on a full grid
//   REF-DFS   Analyzer with REFINEMENT DFS gives the same labels, cluster count, sizes, seeds,
//             seed hits, refined voxels and probes as REFINEMENT SKIP_DFS, for SKIP 1-3,
//             CONNECTIVITY 6/18/26, PBC none and xyz
//   STAGE     labelGrid reports clear/probe/refine times >= 0 whose sum does not exceed elapsed,
//             and names its refinement

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "audit_common.hpp"
#include "bls/BLS.hpp"
#include "bls/Options.hpp"
#include "grid/Grid.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"

namespace {

using Key = std::tuple<int, int, int>;

struct Combo {
  const char* name;
  bls::LatticeType t;
  bls::CenteringType c;
};
const Combo kCombos[] = {
    {"cub-P", bls::LatticeType::Cubic, bls::CenteringType::P},
    {"cub-I", bls::LatticeType::Cubic, bls::CenteringType::I},
    {"cub-F", bls::LatticeType::Cubic, bls::CenteringType::F},
    {"hex-P", bls::LatticeType::Hexagonal, bls::CenteringType::P},
    {"hex-F", bls::LatticeType::Hexagonal, bls::CenteringType::F},
    {"tri-I", bls::LatticeType::Triclinic, bls::CenteringType::I},
};

std::vector<Key> seedList(const bls::Enumerator& en) {
  std::vector<Key> v;
  en.forEach([&](const bls::Enumerator::Seed& s) { v.emplace_back(s.x, s.y, s.z); });
  return v;
}

bls::LatticeSettings settingsFor(bls::LatticeType t, bls::CenteringType c) {
  bls::LatticeSettings s;
  s.lattice = t;
  s.centering = c;
  s.latticeSet = s.centeringSet = true;
  return s;
}

std::string vecStr(const bls::Vec3& v) {
  return "(" + std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z) + ")";
}

// ---------------------------------------------------------------------------
void separableIdentity() {
  std::mt19937_64 rng(280926);
  const double dnns[] = {0.6, 0.9, 1.41, 2.4, 3.3};
  const double ps[] = {0.02, 0.3, 1.0};
  const bls::Vec3 origins[] = {bls::Vec3{}, bls::Vec3{0.3, -0.7, 1.25}, bls::Vec3{1, 2, 3},
                               bls::Vec3{0.5, 0.5, 0.5}, bls::Vec3{-2.5, 7.5, -0.5},
                               bls::Vec3{0.1234567, 0.7654321, 0.4999999999}};
  std::uniform_int_distribution<int> dim(1, 29);
  int open = 0, periodic = 0;
  for (const auto& c : kCombos) {
    const auto lat = settingsFor(c.t, c.c);
    const auto d = bls::buildLattice(lat);
    for (double dnn : dnns)
      for (double p : ps)
        for (const auto& o : origins) {
          const int nx = dim(rng), ny = dim(rng), nz = dim(rng);
          const auto g = audit::randomGrid(rng, nx, ny, nz, p);
          const bls::Mat3 sb = d.basis * (dnn / d.dmin);
          const bls::Enumerator fast(sb, d.offsets, nx, ny, nz, g.occ, o);
          const auto ref = bls::Enumerator::rowWise(sb, d.offsets, nx, ny, nz, g.occ, o);
          ++open;
          const std::string tag = std::string(c.name) + " dNN=" + std::to_string(dnn) + " p=" +
                                  std::to_string(p) + " origin=" + vecStr(o) + " grid " +
                                  std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz);
          AUDIT_CHECK("SEP-ID", seedList(fast) == seedList(ref) && fast.probes() == ref.probes(),
                      tag + ": seeds " + std::to_string(fast.count()) + " vs " + std::to_string(ref.count()) +
                          ", probes " + std::to_string(fast.probes()) + " vs " + std::to_string(ref.probes()));
          if (p == 1.0) {
            const bls::Enumerator all(sb, d.offsets, nx, ny, nz);
            if (o.x == 0.0 && o.y == 0.0 && o.z == 0.0) {
              AUDIT_CHECK("SEP-ALL", seedList(all) == seedList(fast) && all.probes() == fast.probes(),
                          tag + ": lattice-only constructor differs from the full-grid probe set");
            }
          }
          if (c.t != bls::LatticeType::Cubic) continue;
          bls::Mat3 B;
          double eff = 0.0;
          std::string err;
          if (!bls::commensurateCubicBasis(lat, d, dnn, nx, ny, nz, B, eff, err)) continue;
          const bls::Enumerator pf(B, d.offsets, nx, ny, nz, g.occ, bls::Enumerator::PeriodicTag{}, o);
          const auto pr = bls::Enumerator::rowWisePeriodic(B, d.offsets, nx, ny, nz, g.occ, o);
          ++periodic;
          AUDIT_CHECK("SEP-ID", seedList(pf) == seedList(pr) && pf.probes() == pr.probes(),
                      tag + " (PBC): seeds " + std::to_string(pf.count()) + " vs " + std::to_string(pr.count()) +
                          ", probes " + std::to_string(pf.probes()) + " vs " + std::to_string(pr.probes()));
        }
  }
  // Production-sized grid, working deck geometry (cubic F, dNN 1.41 voxels), sparse occupancy.
  {
    const auto lat = settingsFor(bls::LatticeType::Cubic, bls::CenteringType::F);
    const auto d = bls::buildLattice(lat);
    const auto g = audit::randomGrid(rng, 120, 97, 133, 0.004);
    const bls::Mat3 sb = d.basis * (1.41 / d.dmin);
    for (const auto& o : origins) {
      const bls::Enumerator fast(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
      const auto ref = bls::Enumerator::rowWise(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
      ++open;
      AUDIT_CHECK("SEP-ID", seedList(fast) == seedList(ref) && fast.probes() == ref.probes(),
                  "120x97x133 cub-F dNN 1.41 origin " + vecStr(o));
    }
  }
  std::printf("measure: SEP-ID %d open + %d periodic random cases\n", open, periodic);
}

// ---------------------------------------------------------------------------
struct BlsRun {
  bool ok{false};
  std::string err;
  bls::FrameMetrics m;
  std::vector<int> labels;
};
BlsRun runAnalyzer(const audit::Grid3& g, const bls::LatticeSettings& lat, double dnn, int stride,
                   int connectivity, bool periodic, bls::RefinementMode mode, bool withLabels = true) {
  bls::BLSConfig cfg;
  cfg.lattice = lat;
  cfg.hasExplicitDnn = true;
  cfg.dnn = dnn;
  cfg.connectivity = connectivity;
  cfg.refinementStride = stride;
  cfg.refinement = mode;
  cfg.gridSpacing = 1.0;
  bls::Grid grid;
  const bls::Mat3 box{bls::Vec3{double(g.nx), 0, 0}, bls::Vec3{0, double(g.ny), 0},
                      bls::Vec3{0, 0, double(g.nz)}};
  grid.configure(g.nx, g.ny, g.nz, 1.0, box, bls::Vec3{},
                 periodic ? bls::BoxPeriodicity::Periodic : bls::BoxPeriodicity::NonPeriodic);
  grid.occupancy() = g.occ;
  BlsRun r;
  bls::Analyzer an(cfg);
  r.ok = an.labelGrid(grid, r.m, r.err, withLabels ? &r.labels : nullptr);
  return r;
}

void refinementEquivalence() {
  std::mt19937_64 rng(2809261);
  std::uniform_int_distribution<int> dim(3, 26);
  int cases = 0;
  for (const auto& c : {kCombos[0], kCombos[2]})
    for (double p : {0.1, 0.3, 0.5})
      for (double dnn : {0.9, 1.41, 2.4})
        for (int stride = 1; stride <= 3; ++stride)
          for (int conn : {6, 18, 26})
            for (bool periodic : {false, true}) {
              const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
              const auto lat = settingsFor(c.t, c.c);
              const auto a = runAnalyzer(g, lat, dnn, stride, conn, periodic, bls::RefinementMode::SkipDFS);
              const auto b = runAnalyzer(g, lat, dnn, stride, conn, periodic, bls::RefinementMode::DFS);
              ++cases;
              const bool same = a.ok && b.ok && a.labels == b.labels && a.m.nclusters == b.m.nclusters &&
                                a.m.clusterSizes == b.m.clusterSizes && a.m.seeds == b.m.seeds &&
                                a.m.seedHits == b.m.seedHits && a.m.refinedVoxels == b.m.refinedVoxels &&
                                a.m.probes == b.m.probes && a.m.maxCluster == b.m.maxCluster;
              AUDIT_CHECK("REF-DFS", same,
                          std::string(c.name) + " p=" + std::to_string(p) + " dNN=" + std::to_string(dnn) +
                              " SKIP " + std::to_string(stride) + " conn " + std::to_string(conn) +
                              (periodic ? " PBC xyz" : " PBC none") + ": skip_dfs " +
                              std::to_string(a.m.nclusters) + " clusters vs dfs " + std::to_string(b.m.nclusters) +
                              " " + a.err + b.err);
              AUDIT_CHECK("REF-DFS", a.m.refinement == "skip_dfs" && b.m.refinement == "dfs",
                          "refinement names: " + a.m.refinement + " / " + b.m.refinement);
            }
  std::printf("measure: REF-DFS %d cases\n", cases);
}

void stageTimes() {
  std::mt19937_64 rng(2809262);
  const auto g = audit::randomGrid(rng, 60, 60, 60, 0.01);
  const auto lat = settingsFor(bls::LatticeType::Cubic, bls::CenteringType::F);
  for (auto mode : {bls::RefinementMode::SkipDFS, bls::RefinementMode::DFS})
    for (bool periodic : {false, true}) {
      const auto r = runAnalyzer(g, lat, 1.41, 1, 6, periodic, mode, /*withLabels=*/false);
      const auto& m = r.m;
      const bool good = r.ok && m.clearMs >= 0.0 && m.probeMs >= 0.0 && m.refineMs >= 0.0 &&
                        m.clearMs + m.probeMs + m.refineMs <= m.elapsedMs + 1e-9 && !m.refinement.empty();
      AUDIT_CHECK("STAGE", good,
                  "clear " + std::to_string(m.clearMs) + " probe " + std::to_string(m.probeMs) + " refine " +
                      std::to_string(m.refineMs) + " elapsed " + std::to_string(m.elapsedMs) + " " + m.refinement);
    }
}

}  // namespace

int main() {
  separableIdentity();
  refinementEquivalence();
  stageTimes();
  return audit::finish("test_audit_bls_stages");
}
