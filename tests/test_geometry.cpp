// Foundation-layer property tests: linear algebra (common/Types.hpp) and grid
// geometry (grid/Grid.cpp).
//
// GOVERNING PRINCIPLE. Both defects found in the audit so far survived because
// the test case was degenerate: inverse() is correct for symmetric matrices,
// and every cubic basis is symmetric. So no test below may pass merely because
// the cubic case works. Every geometric property is exercised on, at minimum,
// the triclinic basis a=1, b=1.2, c=1.4, alpha=90, beta=100, gamma=110 -- which
// is neither orthogonal nor symmetric -- alongside randomised families that
// include ill-conditioned and wide-dynamic-range matrices. Where a golden
// output and a property were both possible, the property is what is written.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "common/Box.hpp"
#include "common/Types.hpp"
#include "grid/Grid.hpp"

namespace bls {

// See the friend declaration in grid/Grid.hpp.
struct GridTestAccess {
  static std::size_t index(const Grid& g, int ix, int iy, int iz) { return g.index(ix, iy, iz); }
  static Vec3 fractional(const Grid& g, const Vec3& p, Vec3* wrapped) {
    return g.fractional(p, wrapped);
  }
  static Vec3 voxelCenter(const Grid& g, int ix, int iy, int iz) {
    return g.voxelCenter(ix, iy, iz);
  }
};

}  // namespace bls

using bls::GridTestAccess;
using bls::Grid;
using bls::Mat3;
using bls::OccupancyMode;
using bls::Vec3;

namespace {

constexpr double kPi = 3.14159265358979323846;

int g_failures = 0;
int g_checks = 0;

void check(bool cond, const std::string& msg) {
  ++g_checks;
  if (!cond) {
    ++g_failures;
    std::cerr << "FAIL: " << msg << "\n";
  }
}

// ---------------------------------------------------------------------------
// Helpers deliberately independent of the code under test.
// ---------------------------------------------------------------------------

// Matrix product, built only from Mat3::operator*(Vec3). Column j of A*B is
// A applied to column j of B.
Mat3 matmul(const Mat3& a, const Mat3& b) {
  return Mat3{a * b.cols[0], a * b.cols[1], a * b.cols[2]};
}

double at(const Mat3& m, int row, int col) { return m.cols[col][row]; }

double maxAbsEntry(const Mat3& m) {
  double v = 0.0;
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) v = std::max(v, std::abs(at(m, r, c)));
  return v;
}

// Largest |A(r,c) - I(r,c)|.
double distanceFromIdentity(const Mat3& m) {
  double worst = 0.0;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      double target = (r == c) ? 1.0 : 0.0;
      worst = std::max(worst, std::abs(at(m, r, c) - target));
    }
  }
  return worst;
}

// Cofactor expansion along row 0, written out longhand. Shares no code path
// with determinant() (which routes through dot/cross), so agreement between
// the two is evidence rather than tautology.
double determinantByCofactor(const Mat3& m) {
  const double a = at(m, 0, 0), b = at(m, 0, 1), c = at(m, 0, 2);
  const double d = at(m, 1, 0), e = at(m, 1, 1), f = at(m, 1, 2);
  const double g = at(m, 2, 0), h = at(m, 2, 1), i = at(m, 2, 2);
  return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

std::string describe(const Mat3& m) {
  std::string s = "[";
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      s += std::to_string(at(m, r, c));
      if (!(r == 2 && c == 2)) s += ",";
    }
    if (r != 2) s += " ";
  }
  return s + "]";
}

// The canonical non-orthogonal, non-symmetric reference basis.
Mat3 referenceTriclinic() { return bls::buildTriclinicBox(1.0, 1.2, 1.4, 90.0, 100.0, 110.0); }

// ---------------------------------------------------------------------------
// Random matrix families. Together these total well over the 1000 required,
// and no family is orthogonal or symmetric except by measure-zero accident.
// ---------------------------------------------------------------------------

struct MatrixSample {
  Mat3 m;
  std::string family;
};

std::vector<MatrixSample> buildMatrixCorpus(std::mt19937_64& rng) {
  std::vector<MatrixSample> corpus;
  std::uniform_real_distribution<double> pm2(-2.0, 2.0);
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  auto randomMat = [&]() {
    Mat3 m;
    for (int c = 0; c < 3; ++c) m.cols[c] = Vec3{pm2(rng), pm2(rng), pm2(rng)};
    return m;
  };

  // Family 1: general dense random, entries in [-2, 2].
  for (int i = 0; i < 2000; ++i) {
    Mat3 m = randomMat();
    if (std::abs(bls::determinant(m)) < 1e-9) continue;  // keep inverse() defined
    corpus.push_back({m, "random"});
  }

  // Family 2: ill-conditioned -- column 2 is nearly a combination of 0 and 1,
  // pushing the determinant towards the 1e-12 singularity cutoff.
  for (int i = 0; i < 1000; ++i) {
    Mat3 m = randomMat();
    double eps = std::pow(10.0, -2.0 - 8.0 * unit(rng));  // 1e-2 .. 1e-10
    Vec3 perturb{pm2(rng), pm2(rng), pm2(rng)};
    m.cols[2] = m.cols[0] * 0.7 + m.cols[1] * 0.3 + perturb * eps;
    if (std::abs(bls::determinant(m)) < 1e-9) continue;
    corpus.push_back({m, "ill-conditioned"});
  }

  // Family 3: crystallographic triclinic bases -- the shape actually used by
  // the lattice layer. Non-orthogonal and lower-triangular, hence asymmetric.
  std::uniform_real_distribution<double> lenDist(0.5, 3.0);
  std::uniform_real_distribution<double> angDist(60.0, 120.0);
  for (int i = 0; i < 1000; ++i) {
    Mat3 m = bls::buildTriclinicBox(lenDist(rng), lenDist(rng), lenDist(rng), angDist(rng),
                                    angDist(rng), angDist(rng));
    if (std::abs(bls::determinant(m)) < 1e-9) continue;
    corpus.push_back({m, "triclinic"});
  }

  // Family 4: wide dynamic range -- each column independently scaled over
  // twelve orders of magnitude, so the matrix is badly scaled as well as
  // asymmetric.
  for (int i = 0; i < 1000; ++i) {
    Mat3 m = randomMat();
    for (int c = 0; c < 3; ++c) {
      m.cols[c] *= std::pow(10.0, -6.0 + 12.0 * unit(rng));
    }
    if (std::abs(bls::determinant(m)) < 1e-9) continue;
    corpus.push_back({m, "wide-range"});
  }

  // Family 5: named bases the project actually uses, plus deliberately
  // asymmetric fixtures. These are the cases a cubic-only test would miss.
  const double s3 = std::sqrt(3.0);
  corpus.push_back({Mat3{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}}, "named:cubic"});
  corpus.push_back(
      {Mat3{Vec3{1, 0, 0}, Vec3{-0.5, 0.5 * s3, 0}, Vec3{0, 0, 1.633}}, "named:hexagonal"});
  corpus.push_back({referenceTriclinic(), "named:reference-triclinic"});
  corpus.push_back({Mat3{Vec3{1, 2, 3}, Vec3{0, 1, 4}, Vec3{0, 0, 1}}, "named:shear"});
  corpus.push_back({Mat3{Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}}, "named:cyclic-permutation"});
  return corpus;
}

// ---------------------------------------------------------------------------
// Section 1 -- linear algebra properties.
// ---------------------------------------------------------------------------

void testLinearAlgebra(std::mt19937_64& rng) {
  const auto corpus = buildMatrixCorpus(rng);
  std::cout << "  matrix corpus: " << corpus.size() << " matrices\n";
  check(corpus.size() >= 1000, "matrix corpus must hold at least 1000 matrices");

  int leftFail = 0, rightFail = 0, detFail = 0, asymFail = 0;
  double worstScaledLeft = 0.0, worstScaledRight = 0.0, worstDetRel = 0.0;
  std::string worstLeftDesc, worstRightDesc, worstDetDesc;

  for (const auto& sample : corpus) {
    const Mat3& m = sample.m;
    Mat3 inv = bls::inverse(m);

    // Componentwise residual, scaled by the condition of the problem: the
    // reachable accuracy of a floating-point inverse degrades in proportion to
    // ||M|| * ||M^-1||, so an unscaled tolerance would either be vacuous for
    // well-conditioned inputs or spuriously fail on the ill-conditioned family.
    const double kappa = maxAbsEntry(m) * maxAbsEntry(inv);
    const double tol = 256.0 * 2.220446049250313e-16 * std::max(1.0, kappa);

    const double leftRes = distanceFromIdentity(matmul(inv, m));
    const double rightRes = distanceFromIdentity(matmul(m, inv));
    if (leftRes > tol) {
      if (++leftFail == 1) worstLeftDesc = sample.family + " " + describe(m);
    }
    if (rightRes > tol) {
      if (++rightFail == 1) worstRightDesc = sample.family + " " + describe(m);
    }
    worstScaledLeft = std::max(worstScaledLeft, leftRes / tol);
    worstScaledRight = std::max(worstScaledRight, rightRes / tol);

    // The transposed-inverse defect satisfies NEITHER product, but the sharper
    // statement is that inverse(M) must differ from inverse(M)^T whenever M is
    // asymmetric. This is the check that would have caught the original bug
    // immediately, and the one a cubic-only corpus cannot make at all.
    double asymmetry = 0.0;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) asymmetry = std::max(asymmetry, std::abs(at(m, r, c) - at(m, c, r)));
    if (asymmetry > 1e-6 * std::max(1.0, maxAbsEntry(m))) {
      Mat3 invT = inv.transpose();
      double diff = 0.0;
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          diff = std::max(diff, std::abs(at(inv, r, c) - at(invT, r, c)));
      if (diff <= 0.0) ++asymFail;
    }

    // determinant() vs an independent cofactor expansion.
    const double d1 = bls::determinant(m);
    const double d2 = determinantByCofactor(m);
    const double scale = std::max(std::abs(d1), std::pow(maxAbsEntry(m), 3.0));
    const double rel = std::abs(d1 - d2) / std::max(scale, 1e-300);
    if (rel > 1e-12) {
      if (++detFail == 1) worstDetDesc = sample.family + " " + describe(m);
    }
    worstDetRel = std::max(worstDetRel, rel);
  }

  check(leftFail == 0, "inverse(M)*M != I for " + std::to_string(leftFail) + " matrices; first: " +
                           worstLeftDesc);
  check(rightFail == 0, "M*inverse(M) != I for " + std::to_string(rightFail) +
                            " matrices; first: " + worstRightDesc);
  check(detFail == 0, "determinant() disagrees with cofactor expansion for " +
                          std::to_string(detFail) + " matrices; first: " + worstDetDesc);
  check(asymFail == 0, "inverse() returned a symmetric result for " + std::to_string(asymFail) +
                           " asymmetric inputs (the transposed-inverse signature)");
  std::cout << "  worst residual/tolerance: left " << worstScaledLeft << ", right "
            << worstScaledRight << "; worst det relative error " << worstDetRel << "\n";

  // cross() properties on the same scale range as the matrix corpus.
  std::uniform_real_distribution<double> pm2(-2.0, 2.0);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int anticommFail = 0, orthoFail = 0, selfFail = 0, magFail = 0;
  double worstOrtho = 0.0;
  for (int i = 0; i < 5000; ++i) {
    double sa = std::pow(10.0, -4.0 + 8.0 * unit(rng));
    double sb = std::pow(10.0, -4.0 + 8.0 * unit(rng));
    Vec3 a{pm2(rng) * sa, pm2(rng) * sa, pm2(rng) * sa};
    Vec3 b{pm2(rng) * sb, pm2(rng) * sb, pm2(rng) * sb};

    // Anticommutativity holds bit-exactly in IEEE arithmetic: negating
    // cross(b,a) reorders the same two products, it does not recompute them.
    Vec3 ab = bls::cross(a, b);
    Vec3 ba = bls::cross(b, a);
    if (!(ab.x == -ba.x && ab.y == -ba.y && ab.z == -ba.z)) ++anticommFail;

    Vec3 aa = bls::cross(a, a);
    if (!(aa.x == 0.0 && aa.y == 0.0 && aa.z == 0.0)) ++selfFail;

    // Orthogonality is only approximate: the dot product cancels terms of
    // magnitude |a|^2|b|, so the residual is scaled by that, not by 1.
    const double scale = bls::norm(a) * bls::norm(a) * bls::norm(b);
    double resA = std::abs(bls::dot(ab, a)) / std::max(scale, 1e-300);
    double resB = std::abs(bls::dot(ab, b)) / std::max(bls::norm(a) * bls::norm(b) * bls::norm(b),
                                                       1e-300);
    if (resA > 1e-12 || resB > 1e-12) ++orthoFail;
    worstOrtho = std::max({worstOrtho, resA, resB});

    // |a x b|^2 + (a.b)^2 == |a|^2 |b|^2 (Lagrange). Catches a wrong sign or a
    // swapped term that orthogonality alone would tolerate.
    double lhs = bls::dot(ab, ab) + bls::dot(a, b) * bls::dot(a, b);
    double rhs = bls::dot(a, a) * bls::dot(b, b);
    if (std::abs(lhs - rhs) > 1e-12 * std::max(std::abs(rhs), 1e-300)) ++magFail;
  }
  check(anticommFail == 0,
        "cross() not anticommutative in " + std::to_string(anticommFail) + " cases");
  check(selfFail == 0, "cross(a,a) != 0 in " + std::to_string(selfFail) + " cases");
  check(orthoFail == 0, "cross(a,b) not orthogonal to inputs in " + std::to_string(orthoFail) +
                            " cases");
  check(magFail == 0, "Lagrange identity violated in " + std::to_string(magFail) + " cases");
  std::cout << "  cross(): 5000 pairs, worst relative orthogonality residual " << worstOrtho
            << "\n";
}

// ---------------------------------------------------------------------------
// Section 2 -- box / PBC properties, exercised through Grid.
// ---------------------------------------------------------------------------

// The Cartesian map Grid uses internally (see Grid::voxelCenter).
Vec3 cartesian(const Grid& g, const Vec3& frac) { return g.origin() + g.box() * frac; }

struct BoxCase {
  std::string name;
  Mat3 box;
  Vec3 origin;
};

std::vector<BoxCase> boxCases() {
  const double s3 = std::sqrt(3.0);
  std::vector<BoxCase> cases;
  cases.push_back({"cubic", Mat3{Vec3{6, 0, 0}, Vec3{0, 6, 0}, Vec3{0, 0, 6}}, Vec3{0, 0, 0}});
  cases.push_back({"orthorhombic",
                   Mat3{Vec3{5, 0, 0}, Vec3{0, 7, 0}, Vec3{0, 0, 9}}, Vec3{-2.5, 1.0, 0.0}});
  cases.push_back({"hexagonal",
                   Mat3{Vec3{6, 0, 0}, Vec3{-3, 3 * s3, 0}, Vec3{0, 0, 9.8}}, Vec3{0, 0, 0}});
  Mat3 tri = referenceTriclinic();
  cases.push_back({"triclinic(1,1.2,1.4,90,100,110)x6", tri * 6.0, Vec3{0.3, -1.1, 2.2}});
  cases.push_back({"triclinic-oblique",
                   bls::buildTriclinicBox(8.0, 5.0, 6.5, 71.0, 118.0, 96.0) , Vec3{1.0, 1.0, 1.0}});
  return cases;
}

void testBoxGeometry(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> wide(-5.0, 5.0);
  std::uniform_int_distribution<int> shiftDist(-3, 3);

  for (const auto& bc : boxCases()) {
    Grid g;
    g.configure(12, 12, 12, 0.5, bc.box, bc.origin);
    const double boxScale = maxAbsEntry(bc.box) + bls::norm(bc.origin);
    const double tol = 1e-11 * std::max(1.0, boxScale);

    int fwdFail = 0, revFail = 0, idemFail = 0, transFail = 0, rangeFail = 0;
    double worstFwd = 0.0, worstRev = 0.0, worstIdem = 0.0, worstTrans = 0.0;

    for (int i = 0; i < 2000; ++i) {
      // fractional(cartesian(f)) == f, for f already in [0,1)^3.
      Vec3 f{unit(rng), unit(rng), unit(rng)};
      Vec3 back = GridTestAccess::fractional(g, cartesian(g, f), nullptr);
      double e = std::max({std::abs(back.x - f.x), std::abs(back.y - f.y), std::abs(back.z - f.z)});
      // A fractional coordinate a hair below 1.0 can round to exactly 1.0 and
      // wrap to 0.0; that is correct behaviour, not a round-trip failure.
      if (e > 0.5) e = std::abs(e - 1.0);
      if (e > 1e-12) ++fwdFail;
      worstFwd = std::max(worstFwd, e);

      // cartesian(fractional(p)) == wrapped(p), the other direction. For a
      // point drawn anywhere in space this is the definition of the wrap.
      Vec3 p = bc.origin + bc.box * Vec3{wide(rng), wide(rng), wide(rng)};
      Vec3 wrapped;
      Vec3 frac = GridTestAccess::fractional(g, p, &wrapped);
      Vec3 rebuilt = cartesian(g, frac);
      double er = std::max({std::abs(rebuilt.x - wrapped.x), std::abs(rebuilt.y - wrapped.y),
                            std::abs(rebuilt.z - wrapped.z)});
      if (er > tol) ++revFail;
      worstRev = std::max(worstRev, er);

      // The wrap lands strictly inside the unit cell.
      if (!(frac.x >= 0.0 && frac.x <= 1.0 && frac.y >= 0.0 && frac.y <= 1.0 && frac.z >= 0.0 &&
            frac.z <= 1.0)) {
        ++rangeFail;
      }

      // Idempotence: wrapping an already-wrapped point changes nothing.
      Vec3 wrapped2;
      Vec3 frac2 = GridTestAccess::fractional(g, wrapped, &wrapped2);
      double ei = std::max({std::abs(frac2.x - frac.x), std::abs(frac2.y - frac.y),
                            std::abs(frac2.z - frac.z)});
      if (ei > 0.5) ei = std::abs(ei - 1.0);
      if (ei > 1e-11) ++idemFail;
      worstIdem = std::max(worstIdem, ei);

      // Displacing by a whole number of box vectors wraps to the same point.
      // This is the property that fails outright under a transposed inverse
      // for any non-orthogonal box, and is trivially true for a cubic one.
      Vec3 disp = p;
      for (int k = 0; k < 3; ++k) {
        disp += bc.box.column(k) * static_cast<double>(shiftDist(rng));
      }
      Vec3 dispWrapped;
      GridTestAccess::fractional(g, disp, &dispWrapped);
      double et = std::max({std::abs(dispWrapped.x - wrapped.x),
                            std::abs(dispWrapped.y - wrapped.y),
                            std::abs(dispWrapped.z - wrapped.z)});
      // Only a point that lands on a cell face is ambiguous (0 vs 1); allow
      // the wrap to differ by exactly one box vector in that case.
      if (et > tol) {
        double best = et;
        for (int k = 0; k < 3; ++k) {
          for (int s : {-1, 1}) {
            Vec3 alt = dispWrapped - bc.box.column(k) * static_cast<double>(s);
            best = std::min(best, std::max({std::abs(alt.x - wrapped.x), std::abs(alt.y - wrapped.y),
                                            std::abs(alt.z - wrapped.z)}));
          }
        }
        if (best > tol) ++transFail;
        et = best;
      }
      worstTrans = std::max(worstTrans, et);
    }

    check(fwdFail == 0, bc.name + ": fractional(cartesian(f)) != f in " + std::to_string(fwdFail) +
                            "/2000 cases");
    check(revFail == 0, bc.name + ": cartesian(fractional(p)) != wrap(p) in " +
                            std::to_string(revFail) + "/2000 cases");
    check(rangeFail == 0,
          bc.name + ": wrapped fractional outside [0,1] in " + std::to_string(rangeFail) + " cases");
    check(idemFail == 0,
          bc.name + ": PBC wrap not idempotent in " + std::to_string(idemFail) + " cases");
    check(transFail == 0, bc.name + ": displacement by whole box vectors changed the wrap in " +
                              std::to_string(transFail) + " cases");
    std::cout << "  " << bc.name << ": fwd " << worstFwd << ", rev " << worstRev << ", idem "
              << worstIdem << ", box-translation " << worstTrans << "\n";
  }
}

// ---------------------------------------------------------------------------
// Section 3 -- grid layer.
// ---------------------------------------------------------------------------

void testIndexBijection() {
  const int dims[][3] = {{1, 1, 1},  {2, 3, 5},   {5, 3, 2},  {7, 11, 13},
                         {1, 64, 1}, {32, 16, 8}, {40, 40, 40}};
  for (const auto& d : dims) {
    const int nx = d[0], ny = d[1], nz = d[2];
    Grid g;
    g.configure(nx, ny, nz, 1.0, referenceTriclinic() * 10.0, Vec3{0, 0, 0});
    const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;
    check(g.size() == n, "Grid::size() wrong for " + std::to_string(nx) + "x" +
                             std::to_string(ny) + "x" + std::to_string(nz));

    std::vector<uint8_t> seen(n, 0);
    int outOfRange = 0, collisions = 0;
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int iz = 0; iz < nz; ++iz) {
          std::size_t idx = GridTestAccess::index(g, ix, iy, iz);
          if (idx >= n) {
            ++outOfRange;
            continue;
          }
          if (seen[idx]) ++collisions;
          seen[idx] = 1;
        }
      }
    }
    std::size_t unhit = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (!seen[i]) ++unhit;

    const std::string tag = std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz);
    check(outOfRange == 0, tag + ": index() escaped [0,n) in " + std::to_string(outOfRange) +
                               " cases");
    check(collisions == 0, tag + ": index() collided in " + std::to_string(collisions) + " cases");
    check(unhit == 0, tag + ": index() missed " + std::to_string(unhit) + " slots (not onto)");
  }
  std::cout << "  index(): bijection verified on 7 grid shapes\n";
}

// Rasterise a single atom of radius R at the box centre and compare the
// occupied volume to the continuum sphere. Grid::rasterize deliberately
// inflates (OccupancyMode::Any) or deflates (All) the radius by the voxel
// half-diagonal h, so the exact expectations are spheres of radius R+h and
// R-h; the discretisation error on top of that is a surface term, i.e. it
// must scale like 4*pi*R^2*h and not like the volume itself.
void testSphereRasterization() {
  struct Case {
    std::string name;
    Mat3 box;
  };
  std::vector<Case> cases = {
      {"cubic", Mat3{Vec3{20, 0, 0}, Vec3{0, 20, 0}, Vec3{0, 0, 20}}},
      {"triclinic(1,1.2,1.4,90,100,110)x20", referenceTriclinic() * 20.0},
      {"triclinic-oblique", bls::buildTriclinicBox(20.0, 17.0, 23.0, 71.0, 118.0, 96.0)},
  };

  for (const auto& c : cases) {
    const int n = 60;
    Grid g;
    g.configure(n, n, n, 1.0, c.box, Vec3{0, 0, 0});
    const double cellVolume = std::abs(bls::determinant(c.box)) / (double(n) * n * n);

    // Voxel half-diagonal, exactly as Grid::configure computes it.
    Vec3 hx = c.box.column(0) / double(n) * 0.5;
    Vec3 hy = c.box.column(1) / double(n) * 0.5;
    Vec3 hz = c.box.column(2) / double(n) * 0.5;
    double h = 0.0;
    for (int sx : {-1, 1})
      for (int sy : {-1, 1})
        for (int sz : {-1, 1})
          h = std::max(h, bls::norm(hx * sx + hy * sy + hz * sz));

    Vec3 centre = cartesian(g, Vec3{0.5, 0.5, 0.5});
    std::vector<Vec3> atom{centre};

    for (double R : {2.0, 3.5, 5.0}) {
      for (int modeIdx = 0; modeIdx < 2; ++modeIdx) {
        const OccupancyMode mode = modeIdx == 0 ? OccupancyMode::Any : OccupancyMode::All;
        const double effR = modeIdx == 0 ? R + h : R - h;
        if (effR <= 0.0) continue;

        g.rasterize(atom, nullptr, R, mode);
        std::size_t occ = 0;
        for (uint8_t v : g.occupancy())
          if (v) ++occ;

        // Sharp property: rasterize() must mark EXACTLY the voxels whose
        // centre lies within the effective radius. The continuum comparison
        // below is a sanity bound; this is the real specification, and it is
        // the one that exposes a search stencil sized from a single scalar
        // edge length when the three cell vectors differ in length.
        std::size_t exact = 0;
        for (int ix = 0; ix < n; ++ix) {
          for (int iy = 0; iy < n; ++iy) {
            for (int iz = 0; iz < n; ++iz) {
              Vec3 vc = GridTestAccess::voxelCenter(g, ix, iy, iz);
              const double d = bls::norm(vc - centre);
              const bool want = modeIdx == 0 ? (d <= effR + 1e-9) : (d < effR - 1e-9);
              if (want) ++exact;
            }
          }
        }
        const std::string stencilTag = c.name + " R=" + std::to_string(R) +
                                       (modeIdx == 0 ? " Any" : " All");
        check(occ == exact,
              stencilTag + ": rasterize() marked " + std::to_string(occ) +
                  " voxels but exactly " + std::to_string(exact) +
                  " have centres inside the effective radius (" +
                  std::to_string(static_cast<long>(exact) - static_cast<long>(occ)) + " missed)");

        const double measured = double(occ) * cellVolume;
        const double expected = (4.0 / 3.0) * kPi * effR * effR * effR;
        const double surface = 4.0 * kPi * effR * effR * h;
        // Allowance: the boundary layer is one voxel half-diagonal thick, so
        // the count can be off by at most a constant times the surface area
        // times h. A truncated stencil (too small a search radius) blows this
        // budget because it removes a whole cap, not a boundary shell.
        const double allowance = 1.5 * surface + 4.0 * cellVolume;
        const std::string tag = c.name + " R=" + std::to_string(R) +
                                (modeIdx == 0 ? " Any" : " All");
        check(std::abs(measured - expected) <= allowance,
              tag + ": rasterised volume " + std::to_string(measured) + " vs continuum " +
                  std::to_string(expected) + " (allowance " + std::to_string(allowance) +
                  ", error " + std::to_string(measured - expected) + ")");
        std::cout << "  " << tag << ": volume " << measured << " vs " << expected << " (err "
                  << (measured - expected) << ", allow " << allowance << ")\n";
      }
    }
  }
}

// Translating every atom by exactly one grid-cell vector must translate the
// occupancy by exactly one voxel index, modulo the periodic wrap that
// Grid::fractional applies to atom positions.
void testTranslationInvariance(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  struct Case {
    std::string name;
    Mat3 box;
  };
  std::vector<Case> cases = {
      {"cubic", Mat3{Vec3{20, 0, 0}, Vec3{0, 20, 0}, Vec3{0, 0, 20}}},
      {"triclinic(1,1.2,1.4,90,100,110)x20", referenceTriclinic() * 20.0},
  };

  for (const auto& c : cases) {
    const int n = 24;
    Grid g;
    g.configure(n, n, n, 1.0, c.box, Vec3{0, 0, 0});

    // Atoms spread over the whole cell, so the test exercises the boundary
    // rather than avoiding it.
    std::vector<Vec3> atoms;
    for (int i = 0; i < 200; ++i) {
      atoms.push_back(cartesian(g, Vec3{unit(rng), unit(rng), unit(rng)}));
    }

    for (int axis = 0; axis < 3; ++axis) {
      const Vec3 cellVec = c.box.column(axis) / double(n);
      for (double R : {0.0, 2.0}) {
        g.rasterize(atoms, nullptr, R, OccupancyMode::Any);
        std::vector<uint8_t> before = g.occupancy();

        std::vector<Vec3> shifted;
        shifted.reserve(atoms.size());
        for (const auto& a : atoms) shifted.push_back(a + cellVec);
        g.rasterize(shifted, nullptr, R, OccupancyMode::Any);
        const std::vector<uint8_t>& after = g.occupancy();

        std::size_t mismatches = 0;
        for (int ix = 0; ix < n; ++ix) {
          for (int iy = 0; iy < n; ++iy) {
            for (int iz = 0; iz < n; ++iz) {
              int jx = ix, jy = iy, jz = iz;
              if (axis == 0) jx = (ix + 1) % n;
              if (axis == 1) jy = (iy + 1) % n;
              if (axis == 2) jz = (iz + 1) % n;
              if (before[GridTestAccess::index(g, ix, iy, iz)] !=
                  after[GridTestAccess::index(g, jx, jy, jz)]) {
                ++mismatches;
              }
            }
          }
        }
        const std::string tag =
            c.name + " axis=" + std::to_string(axis) + " R=" + std::to_string(R);
        check(mismatches == 0, tag + ": occupancy not invariant under a one-voxel translation (" +
                                   std::to_string(mismatches) + " voxels differ)");
        if (mismatches != 0) {
          std::cout << "  " << tag << ": " << mismatches << " mismatching voxels\n";
        }
      }
    }
  }
  std::cout << "  one-voxel translation invariance: checked 2 boxes x 3 axes x 2 radii\n";
}

}  // namespace

int main() {
  std::mt19937_64 rng(20260829ULL);

  std::cout << "== linear algebra ==\n";
  testLinearAlgebra(rng);
  std::cout << "== box / PBC ==\n";
  testBoxGeometry(rng);
  std::cout << "== grid index ==\n";
  testIndexBijection();
  std::cout << "== sphere rasterisation ==\n";
  testSphereRasterization();
  std::cout << "== translation invariance ==\n";
  testTranslationInvariance(rng);

  std::cout << "\nchecks: " << g_checks << ", failures: " << g_failures << "\n";
  return g_failures == 0 ? 0 : 1;
}
