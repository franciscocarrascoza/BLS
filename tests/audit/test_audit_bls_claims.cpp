// Created: 2026-09-26T12:52+02:00 | by: CC audit (manager) | purpose: BLS-specific manuscript claims — subset, lattice membership, inverse, discretised covering radius, guarantee over a full offset set, probe count, determinism
//
// Guards claims-ledger IDs:
//   C-1/K-3  every BLS component equals exactly one DFS component (no truncation, split or merge)   §4.2 l.670-678
//   M-6      probe set == analytic Bravais lattice sites rounded to voxels (Defect 1 regression)     §2.1
//   M-6T     same for a non-orthogonal, non-symmetric triclinic basis (Defect 2 regression)
//   INV      bls::inverse() is the true inverse of a non-symmetric matrix (Defect 2)
//   M-2      covering radius: measured on the DISCRETISED probe set vs rho = dNN/sqrt2 (cubic F)      §2.1 eq.(covering_dnn)
//   M-2G     guarantee: a digital ball of diameter sqrt2*dNN always contains a probe voxel,
//            worst case over a full offset set (ball centres on a 1/4-voxel mesh spanning one
//            conventional cell, which contains every integer voxel offset)                        §2.1 l.210-212
//   M-7      probe count ~ 4 M^3 / a^3 (cubic F)                                                   §2.4 eq.(probe_count), §4.1 l.612-615
//   R-2      determinism: identical labels on repeated runs                                        §4 l.556-559
//
// Measurements (covering radius, guarantee failure fraction, probe density) are printed
// with "measure:" so they can be quoted; the assertions state the claim as written.
// Working deck: dNN = ALPHA*l0/h = 2.45*2.4/4.0 = 1.47 voxels. Size-floor matrix: dNN = 2.4 voxels.

#include <cmath>
#include <cstdlib>
#include <set>
#include <string>
#include <tuple>

#include "audit_common.hpp"

namespace {

using audit::Grid3;
using Key = std::tuple<int, int, int>;

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

// Probe voxels of a full-occupancy n^3 grid, as the library enumerates them.
std::set<Key> libraryProbes(const bls::LatticeDescriptor& d, double dnn, int n) {
  Grid3 full = audit::emptyGrid(n, n, n);
  std::fill(full.occ.begin(), full.occ.end(), 1);
  bls::Enumerator en(d.basis * (dnn / d.dmin), d.offsets, n, n, n, full.occ);
  std::set<Key> s;
  en.forEach([&](const bls::Enumerator::Seed& p) { s.emplace(p.x, p.y, p.z); });
  return s;
}

// Independent brute force: every site B*(idx + offset) with B = (conventional basis),
// rounded with llround, kept if inside [0,n)^3. idx range generous.
std::set<Key> analyticProbes(const bls::Vec3 a1, const bls::Vec3 a2, const bls::Vec3 a3,
                             const std::vector<bls::Vec3>& offs, int n, int range) {
  std::set<Key> s;
  for (int i = -range; i <= range; ++i)
    for (int j = -range; j <= range; ++j)
      for (int k = -range; k <= range; ++k)
        for (const auto& o : offs) {
          const double fi = i + o.x, fj = j + o.y, fk = k + o.z;
          const double x = a1.x * fi + a2.x * fj + a3.x * fk;
          const double y = a1.y * fi + a2.y * fj + a3.y * fk;
          const double z = a1.z * fi + a2.z * fj + a3.z * fk;
          const long long vx = std::llround(x), vy = std::llround(y), vz = std::llround(z);
          if (vx >= 0 && vy >= 0 && vz >= 0 && vx < n && vy < n && vz < n)
            s.emplace(static_cast<int>(vx), static_cast<int>(vy), static_cast<int>(vz));
        }
  return s;
}

void latticeMembership() {
  const int n = 24;
  for (double dnn : {1.47, 2.4, 3.3}) {
    // Cubic, from the textbook conventional cells (independent of Basis.cpp):
    // P: a = dNN; I: dNN = a*sqrt3/2; F: dNN = a/sqrt2.
    const std::vector<bls::Vec3> P{{0, 0, 0}}, I{{0, 0, 0}, {0.5, 0.5, 0.5}},
        F{{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    struct C { const char* name; bls::CenteringType c; double a; const std::vector<bls::Vec3>* offs; };
    const C cases[] = {{"cubic-P", bls::CenteringType::P, dnn, &P},
                       {"cubic-I", bls::CenteringType::I, 2.0 * dnn / std::sqrt(3.0), &I},
                       {"cubic-F", bls::CenteringType::F, std::sqrt(2.0) * dnn, &F}};
    for (const auto& c : cases) {
      const auto lib = libraryProbes(audit::latticeFor(bls::LatticeType::Cubic, c.c), dnn, n);
      const auto ana = analyticProbes({c.a, 0, 0}, {0, c.a, 0}, {0, 0, c.a}, *c.offs, n,
                                      static_cast<int>(n / c.a) + 3);
      AUDIT_CHECK("M-6", lib == ana, std::string(c.name) + " dNN=" + std::to_string(dnn) +
                  ": probe set != analytic lattice (lib " + std::to_string(lib.size()) + ", analytic " +
                  std::to_string(ana.size()) + ")");
    }
    // Triclinic reference basis a=1,b=1.2,c=1.4, alpha=90, beta=100, gamma=110 (Box.hpp convention),
    // scaled so that the nearest-neighbour distance equals dNN. Brute-force the scale here.
    for (auto cen : {bls::CenteringType::P, bls::CenteringType::I, bls::CenteringType::F}) {
      const bls::LatticeDescriptor d = audit::latticeFor(bls::LatticeType::Triclinic, cen);
      const bls::Mat3 B = d.basis * (dnn / d.dmin);
      // independent check of dmin: nearest distance over all lattice vectors in a 3-shell
      double dmin = 1e9;
      for (int i = -3; i <= 3; ++i) for (int j = -3; j <= 3; ++j) for (int k = -3; k <= 3; ++k)
        for (const auto& o : d.offsets) {
          const bls::Vec3 v = B * bls::Vec3{i + o.x, j + o.y, k + o.z};
          const double l = bls::norm(v);
          if (l > 1e-9) dmin = std::min(dmin, l);
        }
      AUDIT_CHECK("M-6T", std::fabs(dmin - dnn) < 1e-9 * dnn,
                  "triclinic nearest-neighbour distance " + std::to_string(dmin) + " != dNN " + std::to_string(dnn));
      const auto lib = libraryProbes(d, dnn, n);
      const auto ana = analyticProbes(B.column(0), B.column(1), B.column(2), d.offsets, n, 4 * n);
      AUDIT_CHECK("M-6T", lib == ana, "triclinic dNN=" + std::to_string(dnn) + ": probe set != brute-force lattice (lib " +
                  std::to_string(lib.size()) + ", brute " + std::to_string(ana.size()) + ")");
    }
  }
}

void inverseCheck() {
  const bls::Mat3 M{bls::Vec3{1.0, 0.3, -0.2}, bls::Vec3{0.7, 1.2, 0.4}, bls::Vec3{-0.5, 0.9, 1.4}};  // columns, non-symmetric
  const bls::Mat3 Mi = bls::inverse(M);
  double err = 0.0;
  for (int c = 0; c < 3; ++c) {
    const bls::Vec3 e = Mi * M.column(c);  // must be unit vector e_c
    err = std::max({err, std::fabs(e.x - (c == 0)), std::fabs(e.y - (c == 1)), std::fabs(e.z - (c == 2))});
  }
  AUDIT_CHECK("INV", err < 1e-12, "inverse(M)*M != I, max error " + std::to_string(err));
}

// Distance from voxel centre v (integer coordinates, as BLS uses them) to the
// nearest probe voxel, maximised over voxels at least `margin` from every face.
double discretisedCoveringRadius(const std::set<Key>& probes, int n, int margin) {
  std::vector<uint8_t> is(static_cast<std::size_t>(n) * n * n, 0);
  for (const auto& p : probes) is[(std::get<0>(p) * n + std::get<1>(p)) * n + std::get<2>(p)] = 1;
  double worst = 0.0;
  const int R = margin;
  for (int x = margin; x < n - margin; ++x)
    for (int y = margin; y < n - margin; ++y)
      for (int z = margin; z < n - margin; ++z) {
        double best = 1e9;
        for (int dx = -R; dx <= R; ++dx) for (int dy = -R; dy <= R; ++dy) for (int dz = -R; dz <= R; ++dz)
          if (is[((x + dx) * n + (y + dy)) * n + (z + dz)])
            best = std::min(best, std::sqrt(double(dx * dx + dy * dy + dz * dz)));
        worst = std::max(worst, best);
      }
  return worst;
}

void coveringAndGuarantee() {
  for (double dnn : {1.47, 2.4}) {
    for (auto cen : {bls::CenteringType::P, bls::CenteringType::I, bls::CenteringType::F}) {
      const char* cname = cen == bls::CenteringType::P ? "P" : cen == bls::CenteringType::I ? "I" : "F";
      const double rhoCoef = cen == bls::CenteringType::P ? std::sqrt(3.0) / 2.0
                           : cen == bls::CenteringType::I ? std::sqrt(5.0) / (2.0 * std::sqrt(3.0))
                                                          : 1.0 / std::sqrt(2.0);
      const double rho = rhoCoef * dnn;
      const int n = 32, margin = static_cast<int>(std::ceil(rho + 2.0));
      const auto probes = libraryProbes(audit::latticeFor(bls::LatticeType::Cubic, cen), dnn, n);
      const double measured = discretisedCoveringRadius(probes, n, margin);
      std::printf("measure: cubic-%s dNN=%.2f vox: continuous rho=%.4f, discretised covering radius=%.4f (excess %.4f vox)\n",
                  cname, dnn, rho, measured, measured - rho);
      // Provable bound: rounding moves a site by at most sqrt3/2, so the discretised radius <= rho + sqrt3/2.
      AUDIT_CHECK("M-2", measured <= rho + std::sqrt(3.0) / 2.0 + 1e-9,
                  "discretised covering radius exceeds rho + sqrt3/2");

      // Guarantee (claim as written is for cubic F): digital balls of diameter D = 2*rho.
      // Integer voxel offsets of the lattice only translate the rounded probe pattern
      // (llround(s + o) = llround(s) + o), but the conventional edge a is not an integer
      // number of voxels, so the ROUNDED pattern differs from cell to cell. The worst case
      // is therefore taken over ball centres on a 1/4-voxel mesh spanning the whole grid
      // interior, which contains every rounding phase present in the grid.
      std::vector<uint8_t> is(static_cast<std::size_t>(n) * n * n, 0);
      for (const auto& p : probes) is[(std::get<0>(p) * n + std::get<1>(p)) * n + std::get<2>(p)] = 1;
      auto missRate = [&](double D, double step) {
        long tot = 0, miss = 0;
        const double r = D / 2.0;
        const double lo = std::ceil(r) + 2.0, hi = n - std::ceil(r) - 3.0;
        for (double X = lo; X < hi; X += step) for (double Y = lo; Y < hi; Y += step) for (double Z = lo; Z < hi; Z += step) {
          bool hit = false;
          for (int x = static_cast<int>(std::floor(X - r)); x <= static_cast<int>(std::ceil(X + r)) && !hit; ++x)
            for (int y = static_cast<int>(std::floor(Y - r)); y <= static_cast<int>(std::ceil(Y + r)) && !hit; ++y)
              for (int z = static_cast<int>(std::floor(Z - r)); z <= static_cast<int>(std::ceil(Z + r)) && !hit; ++z)
                if ((x - X) * (x - X) + (y - Y) * (y - Y) + (z - Z) * (z - Z) <= r * r + 1e-9 &&
                    is[(x * n + y) * n + z]) hit = true;
          ++tot;
          if (!hit) ++miss;
        }
        return std::make_pair(miss, tot);
      };
      const auto m0 = missRate(2.0 * rho, 0.25);
      double dSafe = 2.0 * rho;
      while (missRate(dSafe, 0.5).first > 0 && dSafe < 2.0 * rho + 4.0) dSafe += 0.1;  // coarse search on a 1/2-voxel mesh
      std::printf("measure: cubic-%s dNN=%.2f vox: digital balls of diameter 2rho=%.3f missed at %ld of %ld centres; smallest diameter never missed = %.2f vox\n",
                  cname, dnn, 2.0 * rho, m0.first, m0.second, dSafe);
      if (cen == bls::CenteringType::F)
        AUDIT_CHECK("M-2G", m0.first == 0, "cubic-F dNN=" + std::to_string(dnn) +
                    ": a digital ball of diameter sqrt2*dNN escaped all probes at " + std::to_string(m0.first) + " centres");
    }
  }
}

void probeCount() {
  const int n = 96;
  for (double dnn : {1.47, 2.4}) {
    const double a = std::sqrt(2.0) * dnn;
    const auto probes = libraryProbes(audit::latticeFor(bls::LatticeType::Cubic, bls::CenteringType::F), dnn, n);
    const double predicted = 4.0 * n * n * n / (a * a * a);
    const double ratio = probes.size() / predicted;
    std::printf("measure: cubic-F dNN=%.2f vox: probe voxels %zu, 4M^3/a^3 = %.0f, ratio %.4f, probe voxels per voxel %.4f\n",
                dnn, probes.size(), predicted, ratio, probes.size() / double(n * n * n));
    AUDIT_CHECK("M-7", std::fabs(ratio - 1.0) < 0.03, "probe count deviates > 3% from 4M^3/a^3 at dNN=" + std::to_string(dnn));
  }
}

void subsetAndDeterminism(int nGrids) {
  std::mt19937_64 rng(0xC0FFEE11ULL);
  for (const auto& lc : kLattices) {
    const bls::LatticeDescriptor d = audit::latticeFor(lc.t, lc.c);
    for (double dnn : {1.47, 2.4, 3.5}) {
      for (int t = 0; t < nGrids; ++t) {
        Grid3 g = audit::randomGrid(rng, 8 + t % 25, 8 + (t * 7) % 25, 8 + (t * 11) % 25, 0.05 + 0.3 * (t % 7) / 6.0);
        std::vector<int> ref;
        audit::runLabelled(bls::ClusterAlgorithm::TraditionalDFS, g, ref);
        const audit::BlsRaw b = audit::runBlsRaw(g, d, dnn);
        const std::string v = audit::subsetViolation(b.labels, ref);
        AUDIT_CHECK("C-1/K-3", v.empty(), std::string(lc.name) + " dNN=" + std::to_string(dnn) + ": " + v);
        AUDIT_CHECK("C-1/K-3", b.nclusters == static_cast<int>(audit::sizesOf(b.labels).size()),
                    "BLS nclusters != number of labelled components");
        if (t % 10 == 0) {
          const audit::BlsRaw b2 = audit::runBlsRaw(g, d, dnn);
          AUDIT_CHECK("R-2", b2.labels == b.labels && b2.seeds == b.seeds, "BLS not deterministic");
        }
      }
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  const int nGrids = argc > 1 ? std::atoi(argv[1]) : 120;
  inverseCheck();
  latticeMembership();
  probeCount();
  coveringAndGuarantee();
  subsetAndDeterminism(nGrids);
  return audit::finish("test_audit_bls_claims");
}
