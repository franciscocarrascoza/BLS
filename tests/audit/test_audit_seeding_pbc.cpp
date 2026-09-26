// Created: 2026-09-26T17:45+02:00 | by: CC audit (manager) | purpose: WP2 tests — BLS probe-site seeding (decision D3), commensurate periodic probe lattice (D10), periodic refinement, LATTICE_ORIGIN
//
// Guards claims-ledger / finding / decision IDs:
//   D3-ID     the probe evaluation (every lattice site in the grid rounded and read, §2.1.2) returns the
//             identical seed list, in the identical order, as the pre-D3 occupied-voxel sweep, for every
//             lattice x centering, several dNN, occupancies and lattice origins (so counts, sizes and
//             labels of BLS cannot change)
//   D3-M      probes() = number of lattice sites (with multiplicity) whose rounded voxel is in the grid,
//             by brute force over a generous index range
//   D3-BLS    Analyzer::labelGrid (the shipped BLS) gives labels identical to a reference built on the
//             pre-D3 sweep, reports probes() and honours LATTICE_ORIGIN
//   PBC-LAT   D10: n_i = ceil(N_i / a), a_i = N_i / n_i <= a, effective dNN <= dNN; periodic probe count
//             n_x n_y n_z |offsets|; integer lattice origin = periodic translation of the probe set;
//             non-commensurate basis rejected; hexagonal / triclinic probe lattices rejected under PBC
//   PBC-WRAP  a component crossing a periodic face is one cluster under PBC xyz, two without (x, y, z)
//   PBC-ORC   periodic grids: BLS at full coverage == periodic N6 components (union-find oracle);
//             at production dNN every BLS cluster is exactly one oracle component and the cluster count
//             equals the number of oracle components holding a probe voxel (stride 1 and 3)
//   PBC-SMALL grids with 1-3 voxels per axis under PBC (wrap onto itself) give one cluster, no fault

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "audit_common.hpp"
#include "bls/BLS.hpp"
#include "bls/Options.hpp"
#include "grid/Grid.hpp"

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
    {"hex-I", bls::LatticeType::Hexagonal, bls::CenteringType::I},
    {"hex-F", bls::LatticeType::Hexagonal, bls::CenteringType::F},
    {"tri-P", bls::LatticeType::Triclinic, bls::CenteringType::P},
    {"tri-I", bls::LatticeType::Triclinic, bls::CenteringType::I},
    {"tri-F", bls::LatticeType::Triclinic, bls::CenteringType::F},
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

// The shipped BLS on a hand-made occupancy: Analyzer::labelGrid on a configured Grid.
struct BlsRun {
  bool ok{false};
  std::string err;
  bls::FrameMetrics m;
  std::vector<int> labels;
};
BlsRun runAnalyzer(const audit::Grid3& g, const bls::LatticeSettings& lat, double dnn, int stride,
                   bool periodic, bls::Vec3 origin = bls::Vec3{}) {
  bls::BLSConfig cfg;
  cfg.lattice = lat;
  cfg.hasExplicitDnn = true;
  cfg.dnn = dnn;
  cfg.connectivity = 6;
  cfg.refinementStride = stride;
  cfg.gridSpacing = 1.0;
  cfg.latticeOrigin[0] = origin.x;
  cfg.latticeOrigin[1] = origin.y;
  cfg.latticeOrigin[2] = origin.z;
  bls::Grid grid;
  const bls::Mat3 box{bls::Vec3{double(g.nx), 0, 0}, bls::Vec3{0, double(g.ny), 0},
                      bls::Vec3{0, 0, double(g.nz)}};
  grid.configure(g.nx, g.ny, g.nz, 1.0, box, bls::Vec3{},
                 periodic ? bls::BoxPeriodicity::Periodic : bls::BoxPeriodicity::NonPeriodic);
  grid.occupancy() = g.occ;
  BlsRun r;
  bls::Analyzer an(cfg);
  r.ok = an.labelGrid(grid, r.m, r.err, &r.labels);
  return r;
}

// ---------------------------------------------------------------------------
void seedIdentity() {
  std::mt19937_64 rng(260926);
  const double dnns[] = {0.9, 1.47, 2.4, 3.3};
  const double ps[] = {0.05, 0.3, 0.8};
  const bls::Vec3 origins[] = {bls::Vec3{}, bls::Vec3{0.3, -0.7, 1.25}, bls::Vec3{1, 2, 3}};
  std::uniform_int_distribution<int> dim(9, 23);
  int cases = 0;
  for (const auto& c : kCombos) {
    const auto d = audit::latticeFor(c.t, c.c);
    for (double dnn : dnns)
      for (double p : ps)
        for (const auto& o : origins) {
          const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
          const bls::Mat3 sb = d.basis * (dnn / d.dmin);
          const bls::Enumerator probe(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
          const auto sweep = bls::Enumerator::occupiedSweep(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
          const auto a = seedList(probe), b = seedList(sweep);
          ++cases;
          AUDIT_CHECK("D3-ID", a == b,
                      std::string(c.name) + " dNN=" + std::to_string(dnn) + " p=" + std::to_string(p) +
                          ": probe seeds " + std::to_string(a.size()) + " vs sweep " + std::to_string(b.size()) +
                          (a.size() == b.size() ? " (same count, different order/content)" : ""));
          AUDIT_CHECK("D3-ID", probe.probes() > 0 && sweep.probes() == 0,
                      std::string(c.name) + ": probes() " + std::to_string(probe.probes()) + " / sweep " +
                          std::to_string(sweep.probes()));
        }
  }
  std::printf("measure: D3-ID %d random cases (9 lattices x 4 dNN x 3 occupancies x 3 origins)\n", cases);
}

// Brute force: every site B(idx + offset) + origin with |idx_k| <= range, counted when its llround is in the grid.
long long bruteProbeCount(const bls::Mat3& B, const std::vector<bls::Vec3>& offs, int nx, int ny, int nz,
                          const bls::Vec3& origin, int range) {
  long long n = 0;
  for (const auto& o : offs)
    for (int i = -range; i <= range; ++i)
      for (int j = -range; j <= range; ++j)
        for (int k = -range; k <= range; ++k) {
          const bls::Vec3 s = B * (bls::Vec3{double(i), double(j), double(k)} + o) + origin;
          const long long x = std::llround(s.x), y = std::llround(s.y), z = std::llround(s.z);
          if (x >= 0 && y >= 0 && z >= 0 && x < nx && y < ny && z < nz) ++n;
        }
  return n;
}

void probeCount() {
  for (const auto& c : kCombos) {
    const auto d = audit::latticeFor(c.t, c.c);
    for (double dnn : {1.47, 2.4})
      for (const auto& o : {bls::Vec3{}, bls::Vec3{0.4, -1.3, 0.75}}) {
        const int nx = 12, ny = 10, nz = 14;
        std::vector<uint8_t> full(static_cast<std::size_t>(nx) * ny * nz, 1);
        const bls::Mat3 sb = d.basis * (dnn / d.dmin);
        const bls::Enumerator en(sb, d.offsets, nx, ny, nz, full, o);
        const long long brute = bruteProbeCount(sb, d.offsets, nx, ny, nz, o, 40);
        AUDIT_CHECK("D3-M", en.probes() == brute,
                    std::string(c.name) + " dNN=" + std::to_string(dnn) + ": probes() " +
                        std::to_string(en.probes()) + " != brute force " + std::to_string(brute));
      }
  }
}

void analyzerMatchesSweep() {
  std::mt19937_64 rng(1729);
  std::uniform_int_distribution<int> dim(10, 26);
  for (const auto& c : kCombos) {
    const auto d = audit::latticeFor(c.t, c.c);
    const auto lat = settingsFor(c.t, c.c);
    for (double dnn : {1.47, 2.4})
      for (int stride : {1, 3})
        for (const auto& o : {bls::Vec3{}, bls::Vec3{0.5, 1.5, -0.25}}) {
          const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), 0.3);
          const auto r = runAnalyzer(g, lat, dnn, stride, false, o);
          // Reference: the pre-D3 sweep's seeds, same refinement, labels dense in seed order.
          const bls::Mat3 sb = d.basis * (dnn / d.dmin);
          const auto sweep = bls::Enumerator::occupiedSweep(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
          std::vector<uint8_t> visited(g.size(), 0);
          std::vector<int> ref(g.size(), -1);
          bls::SkipDFSConfig cfg{g.nx, g.ny, g.nz, 6, stride};
          bls::SkipDFS dfs(cfg, g.occ, visited);
          int nc = 0;
          sweep.forEach([&](const bls::Enumerator::Seed& s) {
            const std::size_t idx = g.index(s.x, s.y, s.z);
            if (!g.occ[idx] || visited[idx]) return;
            if (dfs.runFrom(s.x, s.y, s.z, &ref, nc) > 0) ++nc;
          });
          const bls::Enumerator probe(sb, d.offsets, g.nx, g.ny, g.nz, g.occ, o);
          AUDIT_CHECK("D3-BLS", r.ok && r.labels == ref && r.m.nclusters == nc &&
                                    r.m.seeds == sweep.count() && r.m.probes == probe.probes() && r.m.probes > 0,
                      std::string(c.name) + " dNN=" + std::to_string(dnn) + " stride=" + std::to_string(stride) +
                          ": Analyzer labels/counts differ from the pre-D3 reference (nclusters " +
                          std::to_string(r.m.nclusters) + " vs " + std::to_string(nc) + ", seeds " +
                          std::to_string(r.m.seeds) + " vs " + std::to_string(sweep.count()) + ", probes " +
                          std::to_string(r.m.probes) + ") " + r.err);
        }
  }
}

// ---------------------------------------------------------------------------
void commensurateLattice() {
  const bls::CenteringType cens[] = {bls::CenteringType::P, bls::CenteringType::I, bls::CenteringType::F};
  const int Ns[] = {1, 2, 7, 31, 64, 100};
  for (auto cen : cens) {
    const auto lat = settingsFor(bls::LatticeType::Cubic, cen);
    const auto d = bls::buildLattice(lat);
    for (double dnn : {0.7, 1.47, 2.4, 5.0})
      for (int nx : Ns)
        for (int nz : {3, 64}) {
          const int ny = 31;
          bls::Mat3 B;
          double eff = 0.0;
          std::string err;
          const bool ok = bls::commensurateCubicBasis(lat, d, dnn, nx, ny, nz, B, eff, err);
          const double a = dnn / d.dmin;
          bool good = ok;
          const int dims[3] = {nx, ny, nz};
          for (int k = 0; k < 3 && good; ++k) {
            const double ak = B.column(k)[k];
            const double cells = dims[k] / ak;
            const int n = static_cast<int>(std::llround(cells));
            good = std::fabs(cells - n) < 1e-9 && ak <= a * (1 + 1e-12) &&
                   n == std::max(1, static_cast<int>(std::ceil(dims[k] / a - 1e-9))) &&
                   (n == 1 || (n - 1) * a < dims[k] + 1e-9);
          }
          good = good && eff <= dnn * (1 + 1e-12) && eff > 0.0;
          AUDIT_CHECK("PBC-LAT", good,
                      "cubic-" + bls::centeringToString(cen) + " dNN=" + std::to_string(dnn) + " N=" +
                          std::to_string(nx) + "," + std::to_string(ny) + "," + std::to_string(nz) + ": " + err);
          if (!ok) continue;
          // Periodic probe count and integer-origin translation of the probe set on a full grid.
          std::vector<uint8_t> full(static_cast<std::size_t>(nx) * ny * nz, 1);
          const bls::Enumerator e0(B, d.offsets, nx, ny, nz, full, bls::Enumerator::PeriodicTag{});
          long long expect = static_cast<long long>(d.offsets.size());
          for (int k = 0; k < 3; ++k) expect *= std::llround(dims[k] / B.column(k)[k]);
          AUDIT_CHECK("PBC-LAT", e0.probes() == expect,
                      "periodic probes " + std::to_string(e0.probes()) + " != n_x n_y n_z |offsets| " + std::to_string(expect));
          const bls::Vec3 t{2, -3, 5};
          const bls::Enumerator e1(B, d.offsets, nx, ny, nz, full, bls::Enumerator::PeriodicTag{}, t);
          std::set<Key> shifted, got;
          e0.forEach([&](const bls::Enumerator::Seed& s) {
            shifted.emplace(((s.x + 2) % nx + nx) % nx, ((s.y - 3) % ny + ny) % ny, ((s.z + 5) % nz + nz) % nz);
          });
          e1.forEach([&](const bls::Enumerator::Seed& s) { got.emplace(s.x, s.y, s.z); });
          AUDIT_CHECK("PBC-LAT", got == shifted, "integer LATTICE_ORIGIN is not a periodic translation of the probe set");
        }
  }
  // Non-commensurate basis refused by the periodic enumerator.
  {
    const auto d = audit::latticeFor(bls::LatticeType::Cubic, bls::CenteringType::F);
    std::vector<uint8_t> full(20 * 20 * 20, 1);
    bool threw = false;
    try {
      bls::Enumerator e(d.basis * 2.3, d.offsets, 20, 20, 20, full, bls::Enumerator::PeriodicTag{});
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    AUDIT_CHECK("PBC-LAT", threw, "periodic enumerator accepted a basis that does not tile the grid");
  }
  // Hexagonal / triclinic probe lattices are rejected under PBC, by the basis builder and by BLS.
  for (auto t : {bls::LatticeType::Hexagonal, bls::LatticeType::Triclinic}) {
    const auto lat = settingsFor(t, bls::CenteringType::F);
    const auto d = bls::buildLattice(lat);
    bls::Mat3 B;
    double eff = 0.0;
    std::string err;
    AUDIT_CHECK("PBC-LAT", !bls::commensurateCubicBasis(lat, d, 1.47, 20, 20, 20, B, eff, err) && !err.empty(),
                "non-cubic probe lattice accepted under PBC");
    auto g = audit::emptyGrid(12, 12, 12);
    g.set(3, 3, 3);
    const auto r = runAnalyzer(g, lat, 1.47, 1, true);
    AUDIT_CHECK("PBC-LAT", !r.ok && r.err.find("cubic") != std::string::npos,
                "BLS ran a non-cubic probe lattice under PBC (err: '" + r.err + "')");
  }
}

void wrapAcrossFaces() {
  const auto lat = settingsFor(bls::LatticeType::Cubic, bls::CenteringType::F);
  for (int axis = 0; axis < 3; ++axis)
    for (int stride : {1, 3}) {
      int dims[3] = {9, 8, 7};
      dims[axis] = 20;
      auto g = audit::emptyGrid(dims[0], dims[1], dims[2]);
      for (int u : {0, 1, 2, 17, 18, 19}) {
        int p[3] = {4, 4, 3};
        p[axis] = u;
        g.set(p[0], p[1], p[2]);
      }
      const auto pbc = runAnalyzer(g, lat, 0.7, stride, true);
      const auto open = runAnalyzer(g, lat, 0.7, stride, false);
      AUDIT_CHECK("PBC-WRAP", pbc.ok && pbc.m.nclusters == 1 && pbc.m.maxCluster == 6,
                  "axis " + std::to_string(axis) + " stride " + std::to_string(stride) +
                      ": PBC xyz gave " + std::to_string(pbc.m.nclusters) + " clusters (want 1 of 6) " + pbc.err);
      AUDIT_CHECK("PBC-WRAP", open.ok && open.m.nclusters == 2 && open.m.maxCluster == 3,
                  "axis " + std::to_string(axis) + ": PBC none gave " + std::to_string(open.m.nclusters) +
                      " clusters (want 2 of 3)");
      AUDIT_CHECK("PBC-WRAP", pbc.ok && pbc.m.dnnVoxel <= 0.7 + 1e-12 && open.m.dnnVoxel == 0.7,
                  "dNN reporting: PBC " + std::to_string(pbc.m.dnnVoxel) + ", none " + std::to_string(open.m.dnnVoxel));
    }
}

void periodicOracle() {
  std::mt19937_64 rng(4242);
  std::uniform_int_distribution<int> dim(6, 20);
  int cases = 0;
  for (auto cen : {bls::CenteringType::P, bls::CenteringType::I, bls::CenteringType::F}) {
    const auto lat = settingsFor(bls::LatticeType::Cubic, cen);
    const auto d = bls::buildLattice(lat);
    for (double p : {0.1, 0.25, 0.35})
      for (int stride : {1, 3})
        for (int rep = 0; rep < 4; ++rep) {
          const auto g = audit::randomGrid(rng, dim(rng), dim(rng), dim(rng), p);
          const auto oracle = audit::periodicComponents(g);
          // Full coverage: a <= 1 voxel puts a probe site within half a voxel of every voxel.
          const double dnnFull = 0.7 * d.dmin;
          const auto full = runAnalyzer(g, lat, dnnFull, stride, true);
          AUDIT_CHECK("PBC-ORC", full.ok && audit::canonical(full.labels) == audit::canonical(oracle),
                      "cubic-" + bls::centeringToString(cen) + " full coverage != periodic components " + full.err);
          for (double dnn : {1.47, 2.4, 4.0}) {
            const auto r = runAnalyzer(g, lat, dnn, stride, true);
            const std::string sub = r.ok ? audit::subsetViolation(r.labels, oracle) : r.err;
            // Expected count: oracle components holding at least one probe voxel.
            bls::Mat3 B;
            double eff = 0.0;
            std::string err;
            bls::commensurateCubicBasis(lat, d, dnn, g.nx, g.ny, g.nz, B, eff, err);
            const bls::Enumerator en(B, d.offsets, g.nx, g.ny, g.nz, g.occ, bls::Enumerator::PeriodicTag{});
            std::set<int> hit;
            en.forEach([&](const bls::Enumerator::Seed& s) { hit.insert(oracle[g.index(s.x, s.y, s.z)]); });
            ++cases;
            AUDIT_CHECK("PBC-ORC", r.ok && sub.empty() && r.m.nclusters == static_cast<int>(hit.size()) &&
                                       r.m.probes == en.probes(),
                        "cubic-" + bls::centeringToString(cen) + " dNN=" + std::to_string(dnn) + ": " + sub +
                            " nclusters " + std::to_string(r.m.nclusters) + " vs components with a probe " +
                            std::to_string(hit.size()));
          }
        }
  }
  std::printf("measure: PBC-ORC %d periodic cases at production dNN (plus full-coverage cases)\n", cases);
}

void smallPeriodic() {
  const auto lat = settingsFor(bls::LatticeType::Cubic, bls::CenteringType::F);
  for (int nx : {1, 2, 3})
    for (int ny : {1, 2, 3})
      for (int nz : {1, 3})
        for (int stride : {1, 3}) {
          auto g = audit::emptyGrid(nx, ny, nz);
          std::fill(g.occ.begin(), g.occ.end(), 1);
          const auto r = runAnalyzer(g, lat, 1.47, stride, true);
          AUDIT_CHECK("PBC-SMALL", r.ok && r.m.nclusters == 1 && r.m.maxCluster == nx * ny * nz,
                      std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + " stride " +
                          std::to_string(stride) + ": " + std::to_string(r.m.nclusters) + " clusters " + r.err);
        }
}

}  // namespace

int main() {
  seedIdentity();
  probeCount();
  analyzerMatchesSweep();
  commensurateLattice();
  wrapAcrossFaces();
  periodicOracle();
  smallPeriodic();
  return audit::finish("test_audit_seeding_pbc");
}
