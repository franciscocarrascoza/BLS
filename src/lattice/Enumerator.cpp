#include "lattice/Enumerator.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include "common/Types.hpp"

namespace bls {

namespace {

struct Bounds {
  double minX = std::numeric_limits<double>::infinity();
  double maxX = -std::numeric_limits<double>::infinity();
  double minY = std::numeric_limits<double>::infinity();
  double maxY = -std::numeric_limits<double>::infinity();
  double minZ = std::numeric_limits<double>::infinity();
  double maxZ = -std::numeric_limits<double>::infinity();
};

// std::llround (round half away from zero) without the libm call, for the probe loops
// (audit 27.09.26, finding BLS-PERF: three libm calls per probe site were 30% of a BLS run).
// Exact, not an approximation: t = trunc(x) is exact in a double for |x| < 2^52, and so is
// f = x - t (Sterbenz), so comparing f with +-0.5 decides the rounding exactly as llround
// does. Larger magnitudes (never a voxel coordinate) fall back to std::llround.
inline long long roundHalfAway(double x) {
  if (!(std::fabs(x) < 4503599627370496.0)) return std::llround(x);
  const long long t = static_cast<long long>(x);
  const double f = x - static_cast<double>(t);
  return f >= 0.5 ? t + 1 : (f <= -0.5 ? t - 1 : t);
}

bool seedLess(const Enumerator::Seed& a, const Enumerator::Seed& b) {
  if (a.x != b.x) return a.x < b.x;
  if (a.y != b.y) return a.y < b.y;
  return a.z < b.z;
}

bool seedEqual(const Enumerator::Seed& a, const Enumerator::Seed& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

// Every off-diagonal entry exactly zero: every cubic basis, and every periodic basis
// (commensurateCubicBasis). Such a basis is separable (see build()).
bool isDiagonal(const Mat3& b) {
  const Vec3 c0 = b.column(0), c1 = b.column(1), c2 = b.column(2);
  return c0.y == 0.0 && c0.z == 0.0 && c1.x == 0.0 && c1.z == 0.0 && c2.x == 0.0 && c2.y == 0.0;
}

// The seed list in lexicographic (x,y,z) order, de-duplicated. `runStart` holds the start of
// each offset's run. When every run is already in order (`runsSorted`, the separable path) the
// runs are merged, O(n log k) for k offsets; otherwise the list is sorted, O(n log n). Either
// way the result is the same sequence: seeds with equal keys are equal values, so a sorted
// arrangement of a multiset of seeds is unique.
void finishSeeds(std::vector<Enumerator::Seed>& seeds, std::vector<std::size_t> runStart,
                 bool runsSorted) {
  if (runsSorted) {
    runStart.push_back(seeds.size());
    for (std::size_t r = 1; r + 1 < runStart.size(); ++r) {
      std::inplace_merge(seeds.begin(), seeds.begin() + static_cast<std::ptrdiff_t>(runStart[r]),
                         seeds.begin() + static_cast<std::ptrdiff_t>(runStart[r + 1]), seedLess);
    }
  } else {
    std::sort(seeds.begin(), seeds.end(), seedLess);
  }
  seeds.erase(std::unique(seeds.begin(), seeds.end(), seedEqual), seeds.end());
}

// The separable probe loop: the sites of one offset are the Cartesian product X x Y x Z of
// three per-axis voxel lists (see build()), one entry per lattice index. probes counts every
// site, entries repeated or not. A voxel repeats in a list only when the lattice is finer than
// the grid along that axis (a < 1 voxel), and then several sites round onto one voxel: it is
// read once for all of them, since the seed list is a set (the de-duplication below would drop
// the copies anyway). With each list sorted and repeat-free, the product comes out in
// lexicographic order. The inner loop is one table read and one occupancy read per voxel.
void evaluateProduct(std::vector<int>& X, std::vector<int>& Y, std::vector<int>& Z, int ny, int nz,
                     const std::vector<uint8_t>* occupancy, long long& probes,
                     std::vector<Enumerator::Seed>& out) {
  probes += static_cast<long long>(X.size()) * static_cast<long long>(Y.size()) *
            static_cast<long long>(Z.size());
  for (std::vector<int>* t : {&X, &Y, &Z}) {
    std::sort(t->begin(), t->end());
    t->erase(std::unique(t->begin(), t->end()), t->end());
  }
  const int* zs = Z.data();
  const std::size_t zn = Z.size();
  for (const int vx : X) {
    for (const int vy : Y) {
      if (!occupancy) {
        for (std::size_t k = 0; k < zn; ++k) out.push_back(Enumerator::Seed{vx, vy, zs[k]});
        continue;
      }
      const uint8_t* occRow =
          occupancy->data() + (static_cast<std::size_t>(vx) * static_cast<std::size_t>(ny) +
                               static_cast<std::size_t>(vy)) * static_cast<std::size_t>(nz);
      for (std::size_t k = 0; k < zn; ++k) {
        if (occRow[zs[k]]) out.push_back(Enumerator::Seed{vx, vy, zs[k]});
      }
    }
  }
}

}  // namespace

// Occupancy-driven enumeration (occupiedSweep). Same seed set as the probe
// evaluation below, found by inverting the question: for each OCCUPIED voxel, is
// there a lattice site that rounds to it?
//
// This was BLS's seeding from Task 10 until the 2026-09-26 audit. It is not the
// method the manuscript describes (§2.1.2: the m ~ 4M^3/a^3 probe sites are
// evaluated); it sweeps all M^3 voxels, as DFS's raster does (audit S1-3). Audit
// decision D3 put the described probe evaluation back into BLS; this path is kept
// as the oracle that the probe evaluation returns the identical seed list.
//
// Why the seed set is provably identical, not merely equal in testing:
// llround(site) == v implies |site_k - v_k| <= 0.5 for every k. Writing
// idx = invBasis*site - offset and f0 = invBasis*v - offset gives
// idx - f0 = invBasis*(site - v), so
//     |idx_k - f0_k| <= 0.5 * sum_j |invBasis[k][j]|
// which is the per-axis search half-width below. Every lattice index that
// could round into v therefore lies inside the searched box; nothing outside
// it needs testing. The acceptance test itself is byte-for-byte the volume
// sweep's -- site = basis * (latticeIdx + offset), then llround -- so a
// candidate is admitted on exactly the same arithmetic, with no tolerance and
// no re-derived predicate. The bound is a superset, so widening it (the eps
// below) can only add candidates that the exact test then rejects; it can
// never change the answer.
//
// Emission order needs no sort: voxels are visited in x-major order, which is
// the lexicographic (x,y,z) order the volume path sorts into, and each voxel is
// tested once, which is what its dedup achieves.
void Enumerator::buildFromOccupancy(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                                    int ny, int nz, const std::vector<uint8_t>& occupancy,
                                    const Vec3& origin) {
  const Mat3 invBasis = inverse(basis);

  // Half-width of the candidate box, per lattice-index axis. Derived above; the
  // epsilon covers rounding in the bound's own arithmetic only.
  const double eps = 1e-9;
  double half[3];
  for (int k = 0; k < 3; ++k) {
    Vec3 r = invBasis.row(k);
    half[k] = 0.5 * (std::fabs(r.x) + std::fabs(r.y) + std::fabs(r.z)) + eps;
  }

  std::vector<Vec3> offsetsReal;
  offsetsReal.reserve(offsets.size());
  for (const auto& o : offsets) offsetsReal.push_back(basis * o + origin);

  const std::size_t total = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                            static_cast<std::size_t>(nz);
  const std::size_t planeYZ = static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);

  for (std::size_t base = 0; base < total; ) {
    // Occupancy is a few tenths of a percent dense here, so skip runs of empty
    // voxels a machine word at a time rather than a byte at a time.
    if (base + 8 <= total && (base & 7u) == 0) {
      std::uint64_t word;
      std::memcpy(&word, occupancy.data() + base, 8);
      if (word == 0) { base += 8; continue; }
    }
    if (!occupancy[base]) { ++base; continue; }

    const int vx = static_cast<int>(base / planeYZ);
    const std::size_t rem = base - static_cast<std::size_t>(vx) * planeYZ;
    const int vy = static_cast<int>(rem / static_cast<std::size_t>(nz));
    const int vz = static_cast<int>(rem % static_cast<std::size_t>(nz));
    const Vec3 v{static_cast<double>(vx), static_cast<double>(vy), static_cast<double>(vz)};

    bool admitted = false;
    for (std::size_t oi = 0; oi < offsets.size() && !admitted; ++oi) {
      const Vec3 f0 = invBasis * (v - offsetsReal[oi]);
      const int ixLo = static_cast<int>(std::ceil(f0.x - half[0]));
      const int ixHi = static_cast<int>(std::floor(f0.x + half[0]));
      if (ixLo > ixHi) continue;
      const int iyLo = static_cast<int>(std::ceil(f0.y - half[1]));
      const int iyHi = static_cast<int>(std::floor(f0.y + half[1]));
      if (iyLo > iyHi) continue;
      const int izLo = static_cast<int>(std::ceil(f0.z - half[2]));
      const int izHi = static_cast<int>(std::floor(f0.z + half[2]));
      if (izLo > izHi) continue;

      for (int ix = ixLo; ix <= ixHi && !admitted; ++ix) {
        for (int iy = iyLo; iy <= iyHi && !admitted; ++iy) {
          for (int iz = izLo; iz <= izHi; ++iz) {
            Vec3 latticeIdx{static_cast<double>(ix), static_cast<double>(iy),
                            static_cast<double>(iz)};
            Vec3 site = basis * (latticeIdx + offsets[oi]) + origin;
            if (static_cast<int>(std::llround(site.x)) == vx &&
                static_cast<int>(std::llround(site.y)) == vy &&
                static_cast<int>(std::llround(site.z)) == vz) {
              admitted = true;
              break;
            }
          }
        }
      }
    }
    if (admitted) seeds_.push_back(Seed{vx, vy, vz});
    ++base;
  }
}

void Enumerator::build(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz,
                       const std::vector<uint8_t>* occupancy, const Vec3& origin,
                       bool allowSeparable) {
  Mat3 invBasis = inverse(basis);
  probes_ = 0;
  // Separable basis (audit 28.09.26, finding BLS-PERF2): see the branch in the offset loop.
  const bool separable = allowSeparable && isDiagonal(basis);
  std::vector<std::size_t> runStart;
  std::vector<int> tabX, tabY, tabZ;

  // A site s is kept iff llround(s_k) is in [0, n_k), i.e. iff s lies in the
  // half-open box [-0.5, n_k - 0.5). Bound the lattice indices with the corners
  // of THAT box, not of [0, n_k] -- otherwise the low-side half-voxel shell is
  // outside the hull and has to be recovered by padding. With the corners
  // stated exactly, a padding of 1 covers only the floor/ceil of the hull
  // itself and is not load-bearing: widening it further adds no sites.
  std::vector<Vec3> corners;
  corners.reserve(8);
  for (int i = 0; i <= 1; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        corners.emplace_back(i ? static_cast<double>(nx) - 0.5 : -0.5,
                             j ? static_cast<double>(ny) - 0.5 : -0.5,
                             k ? static_cast<double>(nz) - 0.5 : -0.5);
      }
    }
  }

  const double dims[3] = {static_cast<double>(nx), static_cast<double>(ny),
                          static_cast<double>(nz)};
  const Vec3 c3 = basis.column(2);
  const Vec3 c0s = basis.column(0), c1s = basis.column(1);
  // Rows along z whose z coordinate does not depend on the row (every cubic and hexagonal
  // basis): c3 = (0, 0, c) and c0.z = c1.z = 0. The z list of one offset is then shared by all
  // its rows (the argument of the separable branch below, applied to z alone).
  const bool zShared = allowSeparable && c3.x == 0.0 && c3.y == 0.0 && c0s.z == 0.0 && c1s.z == 0.0;
  long long zCount = 0;  // |Z| with repeats (sites per row), before Z is made repeat-free
  std::vector<Seed> provisional;

  for (const auto& offset : offsets) {
    runStart.push_back(provisional.size());
    Vec3 offsetReal = basis * offset + origin;
    Bounds bnd;
    for (const auto& corner : corners) {
      Vec3 idx = invBasis * (corner - offsetReal);
      bnd.minX = std::min(bnd.minX, idx.x);
      bnd.maxX = std::max(bnd.maxX, idx.x);
      bnd.minY = std::min(bnd.minY, idx.y);
      bnd.maxY = std::max(bnd.maxY, idx.y);
      bnd.minZ = std::min(bnd.minZ, idx.z);
      bnd.maxZ = std::max(bnd.maxZ, idx.z);
    }

    int ixMin = static_cast<int>(std::floor(bnd.minX)) - 1;
    int ixMax = static_cast<int>(std::ceil(bnd.maxX)) + 1;
    int iyMin = static_cast<int>(std::floor(bnd.minY)) - 1;
    int iyMax = static_cast<int>(std::ceil(bnd.maxY)) + 1;
    int izMin = static_cast<int>(std::floor(bnd.minZ)) - 1;
    int izMax = static_cast<int>(std::ceil(bnd.maxZ)) + 1;

    if (zShared) {
      // The z list of this offset, by the row loop's expressions with the row at the box corner;
      // its range [lo, hi] depends only on P.z, the same for every row (see the separable branch).
      const double a0 = static_cast<double>(ixMin) + offset.x;
      const double b0 = static_cast<double>(iyMin) + offset.y;
      const Vec3 P = basis * Vec3{a0, b0, offset.z} + origin;
      double t0 = (-0.5 - P.z) / c3.z, t1 = (dims[2] - 0.5 - P.z) / c3.z;
      if (t0 > t1) std::swap(t0, t1);
      const double lo = std::max(static_cast<double>(izMin), std::floor(t0) - 1.0);
      const double hi = std::min(static_cast<double>(izMax), std::ceil(t1) + 1.0);
      const Vec3 pxy0 = c0s * a0 + c1s * b0;
      tabZ.clear();
      for (int iz = static_cast<int>(lo); lo <= hi && iz <= static_cast<int>(hi); ++iz) {
        const double c = static_cast<double>(iz) + offset.z;
        const long long v = roundHalfAway((pxy0.z + c3.z * c) + origin.z);
        if (v >= 0 && v < nz) tabZ.push_back(static_cast<int>(v));
      }
      zCount = static_cast<long long>(tabZ.size());
    }

    if (separable) {
      // Diagonal basis: the general loop below evaluates site = (c0*a + c1*b + c3*c) + origin
      // with a = ix + offset.x, b = iy + offset.y, c = iz + offset.z, and here every
      // off-diagonal product is +-0. Adding a signed zero returns a nonzero double unchanged
      // and a zero as a zero, which rounds to voxel 0 whatever its sign; so the rounded x
      // depends on ix alone, y on iy alone and z on iz alone, bit for bit. Each coordinate is
      // therefore rounded once per index along its own axis -- by the general loop's own
      // expressions, with the other indices at the box corner -- instead of once per site,
      // and the sites of this offset are the Cartesian product of the three in-range lists.
      // Rows the general loop skips as rowOutside are rows whose x or y rounds outside the
      // grid, which the lists already exclude; its z range [lo, hi] depends only on P.z,
      // which is the same for every row here. Same voxels, same multiplicities, same probes().
      const double a0 = static_cast<double>(ixMin) + offset.x;
      const double b0 = static_cast<double>(iyMin) + offset.y;
      tabX.clear();
      tabY.clear();
      for (int ix = ixMin; ix <= ixMax; ++ix) {
        const Vec3 pxy = c0s * (static_cast<double>(ix) + offset.x) + c1s * b0;
        const long long v = roundHalfAway(pxy.x + origin.x);
        if (v >= 0 && v < nx) tabX.push_back(static_cast<int>(v));
      }
      for (int iy = iyMin; iy <= iyMax; ++iy) {
        const Vec3 pxy = c0s * a0 + c1s * (static_cast<double>(iy) + offset.y);
        const long long v = roundHalfAway(pxy.y + origin.y);
        if (v >= 0 && v < ny) tabY.push_back(static_cast<int>(v));
      }
      evaluateProduct(tabX, tabY, tabZ, ny, nz, occupancy, probes_, provisional);  // Z: above
      continue;
    }
    if (zShared) {
      // Shared z list, read repeat-free (repeats round onto one voxel; probes_ still counts
      // every site through zCount). Rows are not a product here (hexagonal: x depends on ix and
      // iy), so the run is not in order and finishSeeds sorts.
      std::sort(tabZ.begin(), tabZ.end());
      tabZ.erase(std::unique(tabZ.begin(), tabZ.end()), tabZ.end());
    }

    for (int ix = ixMin; ix <= ixMax; ++ix) {
      for (int iy = iyMin; iy <= iyMax; ++iy) {
        // The sites of this row are s(iz) = P + iz*c3. Walk only the iz whose site can
        // round into the grid: s_k in [-0.5, n_k - 0.5] on every axis, solved for iz and
        // widened by one each side. A superset of the accepted sites -- the exact
        // llround test below still decides -- so the result cannot change; it only
        // stops the index box's corners from costing evaluations outside the grid.
        const Vec3 P = basis * Vec3{static_cast<double>(ix) + offset.x,
                                    static_cast<double>(iy) + offset.y, offset.z} +
                       origin;
        double lo = izMin, hi = izMax;
        bool rowOutside = false;
        for (int k = 0; k < 3 && !rowOutside; ++k) {
          if (std::fabs(c3[k]) < 1e-12) {
            rowOutside = P[k] < -0.5 - 1e-6 || P[k] > dims[k] - 0.5 + 1e-6;
            continue;
          }
          double t0 = (-0.5 - P[k]) / c3[k], t1 = (dims[k] - 0.5 - P[k]) / c3[k];
          if (t0 > t1) std::swap(t0, t1);
          lo = std::max(lo, std::floor(t0) - 1.0);
          hi = std::min(hi, std::ceil(t1) + 1.0);
        }
        if (rowOutside || lo > hi) continue;

        // site = basis * (idx + offset) + origin, evaluated as Mat3::operator* does,
        // (c0*a + c1*b) + c2*c, with the iz-independent part hoisted: the same operations
        // in the same order, so every site is bit-identical to the unhoisted form.
        const Vec3 c0 = basis.column(0), c1 = basis.column(1);
        const double a = static_cast<double>(ix) + offset.x;
        const double b = static_cast<double>(iy) + offset.y;
        const Vec3 pxy = c0 * a + c1 * b;
        if (c3.x == 0.0 && c3.y == 0.0) {
          // Rows along z (every cubic basis): c3.x*c and c3.y*c are +-0, and adding a zero
          // leaves a value unchanged, so site.x = pxy.x + origin.x and site.y likewise for the
          // whole row. They are rounded once per row; the sites and the seeds are the ones the
          // general loop below produces, bit for bit.
          const int vx = static_cast<int>(roundHalfAway(pxy.x + origin.x));
          const int vy = static_cast<int>(roundHalfAway(pxy.y + origin.y));
          if (vx < 0 || vy < 0 || vx >= nx || vy >= ny) continue;
          const std::size_t rowBase =
              (static_cast<std::size_t>(vx) * static_cast<std::size_t>(ny) +
               static_cast<std::size_t>(vy)) * static_cast<std::size_t>(nz);
          const uint8_t* occRow = occupancy ? occupancy->data() + rowBase : nullptr;
          if (zShared) {
            probes_ += zCount;
            for (const int vz : tabZ) {
              if (occRow && !occRow[vz]) continue;
              provisional.push_back(Seed{vx, vy, vz});
            }
            continue;
          }
          long long rowProbes = 0;  // local: push_back may alias a member counter
          for (int iz = static_cast<int>(lo); iz <= static_cast<int>(hi); ++iz) {
            const double c = static_cast<double>(iz) + offset.z;
            const int vz = static_cast<int>(roundHalfAway((pxy.z + c3.z * c) + origin.z));
            if (vz < 0 || vz >= nz) continue;
            ++rowProbes;
            if (occRow && !occRow[vz]) continue;
            provisional.push_back(Seed{vx, vy, vz});
          }
          probes_ += rowProbes;
          continue;
        }
        for (int iz = static_cast<int>(lo); iz <= static_cast<int>(hi); ++iz) {
          const double c = static_cast<double>(iz) + offset.z;
          const Vec3 site = (pxy + c3 * c) + origin;
          int vx = static_cast<int>(roundHalfAway(site.x));
          int vy = static_cast<int>(roundHalfAway(site.y));
          int vz = static_cast<int>(roundHalfAway(site.z));
          if (vx < 0 || vy < 0 || vz < 0 || vx >= nx || vy >= ny || vz >= nz) continue;
          ++probes_;
          if (occupancy) {
            // Reject here rather than after the sort: on sparse grids this is
            // the difference between materialising every lattice site in the
            // box and materialising only the seeds.
            std::size_t idx = static_cast<std::size_t>(vx) * static_cast<std::size_t>(ny) *
                                  static_cast<std::size_t>(nz) +
                              static_cast<std::size_t>(vy) * static_cast<std::size_t>(nz) +
                              static_cast<std::size_t>(vz);
            if (!(*occupancy)[idx]) continue;
          }
          provisional.push_back(Seed{vx, vy, vz});
        }
      }
    }
  }

  finishSeeds(provisional, runStart, separable);
  seeds_ = std::move(provisional);
}

void Enumerator::buildPeriodic(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                               int ny, int nz, const std::vector<uint8_t>& occupancy,
                               const Vec3& origin, bool allowSeparable) {
  const int dims[3] = {nx, ny, nz};
  int cells[3];
  for (int k = 0; k < 3; ++k) {
    const Vec3 col = basis.column(k);
    const double a = col[k];
    bool diagonal = a > 0.0;
    for (int j = 0; j < 3; ++j) {
      if (j != k && col[j] != 0.0) diagonal = false;
    }
    const double n = diagonal ? static_cast<double>(dims[k]) / a : 0.0;
    cells[k] = static_cast<int>(std::llround(n));
    if (!diagonal || cells[k] < 1 || std::fabs(n - cells[k]) > 1e-9 * n) {
      throw std::invalid_argument(
          "Enumerator (periodic): the basis is not commensurate with the grid; build it "
          "with commensurateCubicBasis()");
    }
  }

  probes_ = 0;
  std::vector<Seed> provisional;
  auto wrap = [](long long v, int n) {
    if (v >= 0 && v < n) return static_cast<int>(v);  // the common case, no division
    long long r = v % n;
    return static_cast<int>(r < 0 ? r + n : r);
  };
  std::vector<std::size_t> runStart;
  std::vector<int> tabX, tabY, tabZ;
  for (const auto& offset : offsets) {
    runStart.push_back(provisional.size());
    const Vec3 c0 = basis.column(0), c1 = basis.column(1), c2 = basis.column(2);
    if (allowSeparable) {
      // The basis is diagonal (checked above), so the sites are separable exactly as in build():
      // each wrapped coordinate is computed once per index along its own axis, by the row loop's
      // expressions below, and the sites are the product of the three lists (audit BLS-PERF2).
      const double a0 = offset.x, b0 = offset.y;  // ix = iy = 0
      tabX.clear();
      tabY.clear();
      tabZ.clear();
      for (int ix = 0; ix < cells[0]; ++ix) {
        const Vec3 pxy = c0 * (static_cast<double>(ix) + offset.x) + c1 * b0;
        tabX.push_back(wrap(roundHalfAway(pxy.x + origin.x), nx));
      }
      for (int iy = 0; iy < cells[1]; ++iy) {
        const Vec3 pxy = c0 * a0 + c1 * (static_cast<double>(iy) + offset.y);
        tabY.push_back(wrap(roundHalfAway(pxy.y + origin.y), ny));
      }
      const Vec3 pxy0 = c0 * a0 + c1 * b0;
      for (int iz = 0; iz < cells[2]; ++iz) {
        tabZ.push_back(
            wrap(roundHalfAway((pxy0.z + c2.z * (static_cast<double>(iz) + offset.z)) + origin.z), nz));
      }
      evaluateProduct(tabX, tabY, tabZ, ny, nz, &occupancy, probes_, provisional);
      continue;
    }
    for (int ix = 0; ix < cells[0]; ++ix) {
      for (int iy = 0; iy < cells[1]; ++iy) {
        // Hoisted as in build(): bit-identical to basis * (idx + offset) + origin.
        const Vec3 pxy = c0 * (static_cast<double>(ix) + offset.x) +
                         c1 * (static_cast<double>(iy) + offset.y);
        // The periodic basis is diagonal (commensurate cubic, D10): c2 = (0, 0, a_z), so x and
        // y are constant along the row, exactly as in build() (adding +-0 changes nothing).
        const int rx = wrap(roundHalfAway(pxy.x + origin.x), nx);
        const int ry = wrap(roundHalfAway(pxy.y + origin.y), ny);
        const std::size_t rowBase =
            (static_cast<std::size_t>(rx) * static_cast<std::size_t>(ny) +
             static_cast<std::size_t>(ry)) * static_cast<std::size_t>(nz);
        const uint8_t* occRow = occupancy.data() + rowBase;
        for (int iz = 0; iz < cells[2]; ++iz) {
          const int vz = wrap(roundHalfAway((pxy.z + c2.z * (static_cast<double>(iz) + offset.z)) + origin.z), nz);
          if (!occRow[vz]) continue;
          provisional.push_back(Seed{rx, ry, vz});
        }
        probes_ += cells[2];  // every site of the row is evaluated once
      }
    }
  }

  finishSeeds(provisional, runStart, allowSeparable);
  seeds_ = std::move(provisional);
}

Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz) {
  build(basis, offsets, nx, ny, nz, nullptr);
}

// Same lattice sites as the constructor above, restricted to those whose voxel
// is occupied. Seed selection is decided entirely by lattice geometry;
// occupancy only filters.
Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz, const std::vector<uint8_t>& occupancy, const Vec3& origin) {
  build(basis, offsets, nx, ny, nz, &occupancy, origin);
}

Enumerator::Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                       int nz, const std::vector<uint8_t>& occupancy, PeriodicTag,
                       const Vec3& origin) {
  buildPeriodic(basis, offsets, nx, ny, nz, occupancy, origin);
}

Enumerator Enumerator::rowWise(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny,
                               int nz, const std::vector<uint8_t>& occupancy, const Vec3& origin) {
  Enumerator e;
  e.build(basis, offsets, nx, ny, nz, &occupancy, origin, /*allowSeparable=*/false);
  return e;
}

Enumerator Enumerator::rowWisePeriodic(const Mat3& basis, const std::vector<Vec3>& offsets, int nx,
                                       int ny, int nz, const std::vector<uint8_t>& occupancy,
                                       const Vec3& origin) {
  Enumerator e;
  e.buildPeriodic(basis, offsets, nx, ny, nz, occupancy, origin, /*allowSeparable=*/false);
  return e;
}

Enumerator Enumerator::occupiedSweep(const Mat3& basis, const std::vector<Vec3>& offsets,
                                     int nx, int ny, int nz,
                                     const std::vector<uint8_t>& occupancy,
                                     const Vec3& origin) {
  Enumerator e;
  e.buildFromOccupancy(basis, offsets, nx, ny, nz, occupancy, origin);
  return e;
}

Enumerator Enumerator::legacyRadiusSelection(const Mat3& basis, const std::vector<Vec3>& offsets,
                                             int nx, int ny, int nz,
                                             const std::vector<uint8_t>& occupancy) {
  return Enumerator(LegacyRadiusTag{}, basis, offsets, nx, ny, nz, occupancy);
}

// DEFECTIVE — see the note on Enumerator::legacyRadiusSelection in Enumerator.hpp.
// Preserved verbatim so pre-2026-08-29 results can be reproduced. The bug is the
// unconditional per-offset radius test below: `< 0.25` is a fixed fractional
// radius of 0.5 applied once per centering offset, so multi-offset centerings
// capture overlapping spheres covering nearly the entire grid rather than the
// lattice sites. Not to be repaired by retuning the constant — a radius fitted
// to cubic P/I/F would still be wrong for hexagonal and triclinic.
Enumerator::Enumerator(LegacyRadiusTag, const Mat3& basis, const std::vector<Vec3>& offsets,
                       int nx, int ny, int nz, const std::vector<uint8_t>& occupancy) {
  // Frozen copy of the pre-2026-08-29 bls::inverse(), which returned the
  // TRANSPOSED inverse. This constructor exists solely to reproduce pre-fix
  // numbers, so it must reproduce the pre-fix arithmetic exactly -- not merely
  // the pre-fix algebra. Calling the corrected inverse() and transposing it
  // back is not good enough: the radius test below is a strict `< 0.25` on a
  // knife edge that voxels land on exactly (at dNN = 2.4 voxels, every voxel
  // 1.2 from a site), so a one-ulp difference in the matrix flips hundreds of
  // seeds. Do not replace this with a call to inverse().
  auto frozenPreFixInverse = [](const Mat3& m) {
    Vec3 c0 = cross(m.cols[1], m.cols[2]);
    Vec3 c1 = cross(m.cols[2], m.cols[0]);
    Vec3 c2 = cross(m.cols[0], m.cols[1]);
    double det = dot(m.cols[0], c0);
    return Mat3{c0 / det, c1 / det, c2 / det};
  };
  Mat3 invBasis = frozenPreFixInverse(basis);

  // Precompute real-space offsets
  std::vector<Vec3> offsetsReal;
  offsetsReal.reserve(offsets.size());
  for (const auto& o : offsets) {
    offsetsReal.push_back(basis * o);
  }

  // Scan the occupancy grid sequentially (cache-friendly).
  // For each occupied voxel, check if it lies on any lattice site.
  for (int x = 0; x < nx; ++x) {
    std::size_t xStride = static_cast<std::size_t>(x) * static_cast<std::size_t>(ny) *
                          static_cast<std::size_t>(nz);
    for (int y = 0; y < ny; ++y) {
      std::size_t xyStride = xStride + static_cast<std::size_t>(y) * static_cast<std::size_t>(nz);
      for (int z = 0; z < nz; ++z) {
        if (!occupancy[xyStride + static_cast<std::size_t>(z)]) continue;

        Vec3 pos{static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)};
        for (const auto& oReal : offsetsReal) {
          Vec3 frac = invBasis * (pos - oReal);
          double dx = frac.x - std::round(frac.x);
          double dy = frac.y - std::round(frac.y);
          double dz = frac.z - std::round(frac.z);
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            seeds_.push_back(Seed{x, y, z});
            break;  // one match is enough
          }
        }
      }
    }
  }
}

}  // namespace bls

