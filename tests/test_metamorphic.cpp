// Task 5: metamorphic invariants, applied to every algorithm rather than just BLS.
//
// A metamorphic test does not need to know the right answer. It transforms an
// input in a way whose effect on the output is known -- translate the pattern
// and the component count must not change; add one voxel and the count can rise
// by at most one -- and checks the relation. That makes it the only practical
// way to test the seven approximate methods here, for which no ground truth
// exists at all.
//
// Every invariant is run against all fourteen methods. Some are EXPECTED to
// fail for the approximate ones: k-means partitions occupied space into k
// pieces regardless of connectivity, so "adding one voxel changes the count by
// at most one" is not a claim it makes. The table therefore separates
// violations that indicate a defect from violations that are the method
// behaving as designed, and the expectation for each pair is declared up front
// rather than read off the results afterwards.
//
// BLS is included at two seed densities. Its lattice is anchored to the grid
// origin, so translating a pattern by a non-lattice vector changes which
// components a seed lands in. Whether that shows up as a translation-invariance
// violation is exactly the coverage question Task 6 measures, and it is worth
// seeing the size of the effect here first.
//
// The "translation under PBC" invariant from the original plan is replaced by
// translation strictly inside a larger grid: with no periodic boundary anywhere
// in this codebase, a wrapping translation is not a transformation any of these
// methods is defined under.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"
#include "refine/SkipDFS.hpp"

using bls::ClusterAlgorithm;
using bls::ClusterParams;
using bls::ClusterResult;

namespace {

struct Grid3 {
  int nx{0}, ny{0}, nz{0};
  std::vector<uint8_t> occ;
  std::size_t index(int i, int j, int k) const {
    return static_cast<std::size_t>(i) * ny * nz + static_cast<std::size_t>(j) * nz + k;
  }
  std::size_t size() const { return static_cast<std::size_t>(nx) * ny * nz; }
  long occupied() const {
    long n = 0;
    for (uint8_t v : occ) if (v) ++n;
    return n;
  }
};

Grid3 makeGrid(int nx, int ny, int nz) {
  Grid3 g; g.nx = nx; g.ny = ny; g.nz = nz; g.occ.assign(g.size(), 0); return g;
}

// --- reference labelling, for invariants that need to know the components ---

std::vector<int> referenceLabels(const Grid3& g) {
  static const int d6[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  std::vector<int> labels(g.size(), -1);
  std::vector<std::size_t> stack;
  int next = 0;
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k) {
        std::size_t s = g.index(i, j, k);
        if (!g.occ[s] || labels[s] >= 0) continue;
        const int id = next++;
        stack.clear(); stack.push_back(s); labels[s] = id;
        while (!stack.empty()) {
          std::size_t cur = stack.back(); stack.pop_back();
          int ci = static_cast<int>(cur / (static_cast<std::size_t>(g.ny) * g.nz));
          int rem = static_cast<int>(cur % (static_cast<std::size_t>(g.ny) * g.nz));
          int cj = rem / g.nz, ck = rem % g.nz;
          for (const auto& d : d6) {
            int ni = ci + d[0], nj = cj + d[1], nk = ck + d[2];
            if (ni < 0 || nj < 0 || nk < 0 || ni >= g.nx || nj >= g.ny || nk >= g.nz) continue;
            std::size_t nidx = g.index(ni, nj, nk);
            if (!g.occ[nidx] || labels[nidx] >= 0) continue;
            labels[nidx] = id; stack.push_back(nidx);
          }
        }
      }
  return labels;
}

// --- the methods under test -------------------------------------------------

struct Outcome {
  int nclusters{0};
  int maxCluster{0};
  long sumSizes{0};
  std::vector<int> sizes;
};

struct Method {
  std::string name;
  ClusterAlgorithm algo{ClusterAlgorithm::BLS};
  bool isBls{false};
  double blsDnnVoxels{0.0};
  bool exact{false};  // claims an exact connected-component partition

  // True when the method's output depends on where the pattern sits relative
  // to the grid origin, or on the order voxels are scanned in -- in which case
  // translation and rotation equivariance are not claims it makes. Set from
  // how each method is implemented, not from how it scored:
  //   bls   - lattice sites are enumerated from the grid origin
  //           (Enumerator walks lattice indices over the grid volume)
  //   vccs  - seeds placed at si = halfS; si < nx; si += S, an origin-anchored
  //           uniform lattice (Algorithms.cpp, "Place seeds on a regular 3D grid")
  //   kmeans- centroids initialised at points[(i*numPoints)/k], indices into a
  //           list built by scanning the grid in (i,j,k) order, so the
  //           initialisation rotates with the scan and not with the geometry
  bool anchored{false};
};

std::vector<Method> methods() {
  auto m = [](ClusterAlgorithm a, bool exact) {
    Method x; x.algo = a; x.name = bls::algorithmToString(a); x.exact = exact; return x;
  };
  std::vector<Method> v = {
      m(ClusterAlgorithm::TraditionalDFS, true),
      m(ClusterAlgorithm::CC3D, true),
      m(ClusterAlgorithm::CC3DOptimized, true),
      m(ClusterAlgorithm::GCBD, true),
      m(ClusterAlgorithm::RLECCL, true),
      m(ClusterAlgorithm::RLECCLOptimized, true),
      m(ClusterAlgorithm::SkipDFS, false),
      m(ClusterAlgorithm::DBSCAN, false),
      m(ClusterAlgorithm::Hierarchical, false),
      m(ClusterAlgorithm::KMeans, false),   // .anchored set below
      m(ClusterAlgorithm::VCCSOptimized, false),
      m(ClusterAlgorithm::HDBSCAN, false),
      m(ClusterAlgorithm::VCCS, false),
  };
  for (auto& x : v) {
    if (x.algo == ClusterAlgorithm::KMeans || x.algo == ClusterAlgorithm::VCCS ||
        x.algo == ClusterAlgorithm::VCCSOptimized) {
      x.anchored = true;
    }
  }
  Method dense; dense.name = "bls(dNN=1.5vox)"; dense.isBls = true; dense.blsDnnVoxels = 1.5; dense.anchored = true;
  Method sparse; sparse.name = "bls(dNN=3.0vox)"; sparse.isBls = true; sparse.blsDnnVoxels = 3.0; sparse.anchored = true;
  v.push_back(dense);
  v.push_back(sparse);
  return v;
}

// BLS at the grid level is exactly what Analyzer does once rasterisation is
// done: enumerate lattice sites filtered by occupancy, then Skip-DFS from each.
// Reproduced here so BLS can be driven with a synthetic occupancy grid, which
// the Analyzer API cannot do (it takes atom positions).
Outcome runBls(const Grid3& g, double dnnVoxels) {
  bls::LatticeSettings settings;
  settings.lattice = bls::LatticeType::Cubic;
  settings.centering = bls::CenteringType::F;
  bls::LatticeDescriptor desc = bls::buildLattice(settings);
  bls::Mat3 scaled = desc.basis * (dnnVoxels / desc.dmin);

  std::vector<uint8_t> visited(g.size(), 0);
  bls::Enumerator enumerator(scaled, desc.offsets, g.nx, g.ny, g.nz, g.occ);
  bls::SkipDFSConfig cfg{g.nx, g.ny, g.nz, 6, 1};
  bls::SkipDFS dfs(cfg, g.occ, visited);

  Outcome o;
  enumerator.forEach([&](const bls::Enumerator::Seed& s) {
    int size = dfs.runFrom(s.x, s.y, s.z);
    if (size > 0) {
      ++o.nclusters;
      o.maxCluster = std::max(o.maxCluster, size);
      o.sizes.push_back(size);
      o.sumSizes += size;
    }
  });
  std::sort(o.sizes.begin(), o.sizes.end(), std::greater<int>());
  return o;
}

Outcome run(const Method& m, const Grid3& g) {
  if (m.isBls) return runBls(g, m.blsDnnVoxels);
  ClusterParams p;
  p.nx = g.nx; p.ny = g.ny; p.nz = g.nz;
  p.connectivity = 6;
  p.skipDfsJumpDistance = 1;
  std::vector<uint8_t> visited(g.size(), 0);
  ClusterResult r = bls::runClusterAlgorithm(m.algo, p, g.occ, visited, nullptr);
  Outcome o;
  o.nclusters = r.nclusters;
  o.maxCluster = r.maxCluster;
  o.sizes = r.clusterSizes;
  std::sort(o.sizes.begin(), o.sizes.end(), std::greater<int>());
  for (int s : o.sizes) o.sumSizes += s;
  return o;
}

// --- grid transformations ---------------------------------------------------

Grid3 embed(const Grid3& src, int pad, int dx, int dy, int dz) {
  Grid3 out = makeGrid(src.nx + 2 * pad, src.ny + 2 * pad, src.nz + 2 * pad);
  for (int i = 0; i < src.nx; ++i)
    for (int j = 0; j < src.ny; ++j)
      for (int k = 0; k < src.nz; ++k)
        if (src.occ[src.index(i, j, k)])
          out.occ[out.index(i + pad + dx, j + pad + dy, k + pad + dz)] = 1;
  return out;
}

// Rotation by 90 degrees about an axis, and mirror reflections. Cubic grids
// only, so the result stays in bounds.
Grid3 rotate90(const Grid3& g, int axis) {
  Grid3 out = makeGrid(g.nx, g.ny, g.nz);
  const int n = g.nx;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        if (!g.occ[g.index(i, j, k)]) continue;
        if (axis == 0) out.occ[out.index(i, n - 1 - k, j)] = 1;
        else if (axis == 1) out.occ[out.index(k, j, n - 1 - i)] = 1;
        else out.occ[out.index(n - 1 - j, i, k)] = 1;
      }
  return out;
}

Grid3 mirror(const Grid3& g, int axis) {
  Grid3 out = makeGrid(g.nx, g.ny, g.nz);
  for (int i = 0; i < g.nx; ++i)
    for (int j = 0; j < g.ny; ++j)
      for (int k = 0; k < g.nz; ++k)
        if (g.occ[g.index(i, j, k)]) {
          int a = axis == 0 ? g.nx - 1 - i : i;
          int b = axis == 1 ? g.ny - 1 - j : j;
          int c = axis == 2 ? g.nz - 1 - k : k;
          out.occ[out.index(a, b, c)] = 1;
        }
  return out;
}

// A and B side by side along x with one empty plane between them.
Grid3 concatWithGap(const Grid3& a, const Grid3& b) {
  Grid3 out = makeGrid(a.nx + 1 + b.nx, std::max(a.ny, b.ny), std::max(a.nz, b.nz));
  for (int i = 0; i < a.nx; ++i)
    for (int j = 0; j < a.ny; ++j)
      for (int k = 0; k < a.nz; ++k)
        if (a.occ[a.index(i, j, k)]) out.occ[out.index(i, j, k)] = 1;
  for (int i = 0; i < b.nx; ++i)
    for (int j = 0; j < b.ny; ++j)
      for (int k = 0; k < b.nz; ++k)
        if (b.occ[b.index(i, j, k)]) out.occ[out.index(a.nx + 1 + i, j, k)] = 1;
  return out;
}

// --- bookkeeping ------------------------------------------------------------

enum Inv {
  kTranslation, kRotation, kReflection, kSumSizes, kMaxSize,
  kAddVoxel, kSeparatedPlane, kRemoveLargest, kInvCount
};

const char* kInvName[kInvCount] = {
    "translation", "rotation90", "reflection", "sum==N_occ", "max<=N_occ",
    "add-voxel", "separated", "remove-largest"};

// Whether the invariant is a claim the method actually makes. An approximate
// method failing an exactness invariant is information, not a defect, and is
// counted separately so the two are never confused.
bool expected(const Method& m, Inv inv) {
  switch (inv) {
    case kMaxSize:
      return true;  // universal: a component cannot exceed the occupied set
    case kTranslation:
    case kRotation:
    case kReflection:
      // A method whose result depends only on the connectivity structure of
      // the occupied set must be equivariant under a symmetry of the voxel
      // lattice. A method that seeds from the grid origin, or initialises from
      // scan order, is not -- by construction, not by defect. See
      // Method::anchored.
      return !m.anchored;
    case kSumSizes:
    case kAddVoxel:
    case kSeparatedPlane:
    case kRemoveLargest:
      return m.exact;
    default:
      return false;
  }
}

struct Cell { long checked = 0; long violated = 0; };

}  // namespace

int main(int argc, char** argv) {
  const int nBase = argc > 1 ? std::stoi(argv[1]) : 60;

  const auto ms = methods();
  std::map<std::string, std::array<Cell, kInvCount>> table;
  for (const auto& m : ms) table[m.name];

  auto record = [&](const Method& m, Inv inv, bool ok) {
    auto& c = table[m.name][inv];
    ++c.checked;
    if (!ok) ++c.violated;
  };

  // Base grids: small, so the O(n^2) approximate methods stay tractable, and
  // cubic, so 90-degree rotations stay in bounds.
  std::vector<Grid3> bases;
  std::mt19937_64 rng(20260829ULL);
  std::uniform_real_distribution<double> u01(0, 1);
  const int sizes[] = {6, 8, 10, 12};
  const double occs[] = {0.08, 0.15, 0.25, 0.3116, 0.40, 0.55};
  for (int t = 0; t < nBase; ++t) {
    const int n = sizes[t % 4];
    const double p = occs[(t / 4) % 6];
    Grid3 g = makeGrid(n, n, n);
    for (auto& v : g.occ) v = u01(rng) < p ? 1 : 0;
    bases.push_back(g);
  }
  // Structured additions, so the sweep is not purely random.
  { Grid3 g = makeGrid(9, 9, 9);
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k)
      g.occ[g.index(i, j, k)] = 1;
    for (int i = 5; i < 9; ++i) for (int j = 5; j < 9; ++j) for (int k = 5; k < 9; ++k)
      g.occ[g.index(i, j, k)] = 1;
    bases.push_back(g); }
  { Grid3 g = makeGrid(9, 9, 9);
    for (int i = 0; i < 9; ++i) g.occ[g.index(i, i, i)] = 1;
    bases.push_back(g); }
  { Grid3 g = makeGrid(8, 8, 8); std::fill(g.occ.begin(), g.occ.end(), 1); bases.push_back(g); }

  std::printf("== metamorphic invariants ==\n");
  std::printf("  %zu base grids x %zu methods\n\n", bases.size(), ms.size());

  std::uniform_int_distribution<int> pick(0, 1 << 30);

  for (std::size_t bi = 0; bi < bases.size(); ++bi) {
    const Grid3& base = bases[bi];
    const long nocc = base.occupied();

    // Precompute the transformed grids once, shared across methods.
    const Grid3 embedA = embed(base, 2, 0, 0, 0);
    const Grid3 embedB = embed(base, 2, 1, -2, 2);
    const Grid3 rot[3] = {rotate90(base, 0), rotate90(base, 1), rotate90(base, 2)};
    const Grid3 mir[3] = {mirror(base, 0), mirror(base, 1), mirror(base, 2)};

    // One extra voxel, chosen among the currently empty ones.
    Grid3 plus = base;
    {
      std::vector<std::size_t> empties;
      for (std::size_t i = 0; i < plus.occ.size(); ++i)
        if (!plus.occ[i]) empties.push_back(i);
      if (!empties.empty()) plus.occ[empties[pick(rng) % empties.size()]] = 1;
    }

    // Two copies separated by an empty plane.
    const Grid3 pair = concatWithGap(base, base);

    // Base with its largest 6-connected component deleted, identified by the
    // reference labelling so the choice does not depend on any method.
    Grid3 minusLargest = base;
    bool hasComponent = false;
    {
      const std::vector<int> lab = referenceLabels(base);
      std::map<int, int> counts;
      for (int v : lab) if (v >= 0) ++counts[v];
      if (!counts.empty()) {
        hasComponent = true;
        int best = counts.begin()->first;
        for (const auto& kv : counts) if (kv.second > counts[best]) best = kv.first;
        for (std::size_t i = 0; i < lab.size(); ++i) if (lab[i] == best) minusLargest.occ[i] = 0;
      }
    }

    for (const auto& m : ms) {
      const Outcome o = run(m, base);

      // 1. Translation inside a larger grid (no PBC exists to wrap under).
      {
        const Outcome a = run(m, embedA);
        const Outcome b = run(m, embedB);
        record(m, kTranslation, a.nclusters == b.nclusters && a.sizes == b.sizes);
      }

      // 2/3. 90-degree rotation and mirror reflection.
      for (int ax = 0; ax < 3; ++ax) {
        const Outcome r = run(m, rot[ax]);
        record(m, kRotation, r.nclusters == o.nclusters && r.sizes == o.sizes);
        const Outcome f = run(m, mir[ax]);
        record(m, kReflection, f.nclusters == o.nclusters && f.sizes == o.sizes);
      }

      // 4. Component sizes must account for every occupied voxel.
      record(m, kSumSizes, o.sumSizes == nocc);

      // 5. No component can be larger than the occupied set.
      record(m, kMaxSize, o.maxCluster <= nocc);

      // 6. Adding one voxel: at most +1 new component; at most 5 destroyed,
      //    since one voxel has 6 neighbours and joining k of them into one
      //    component removes k-1 <= 5 components.
      if (plus.occupied() != nocc) {
        const Outcome p = run(m, plus);
        const int delta = p.nclusters - o.nclusters;
        record(m, kAddVoxel, delta <= 1 && delta >= -5);
      }

      // 7. Two copies separated by an empty plane: counts add.
      {
        const Outcome pr = run(m, pair);
        record(m, kSeparatedPlane, pr.nclusters == 2 * o.nclusters);
      }

      // 8. Removing the largest component drops the count by exactly one.
      if (hasComponent) {
        const Outcome r = run(m, minusLargest);
        record(m, kRemoveLargest, r.nclusters == o.nclusters - 1);
      }
    }
  }

  // --- report ---------------------------------------------------------------
  std::printf("Violations / checks. A dash marks a pair where the invariant is not a claim\n");
  std::printf("the method makes; the count is still shown, in parentheses, as information.\n\n");
  std::printf("  %-18s", "method");
  for (int i = 0; i < kInvCount; ++i) std::printf(" %14s", kInvName[i]);
  std::printf("\n");

  long realViolations = 0;
  for (const auto& m : ms) {
    std::printf("  %-18s", m.name.c_str());
    for (int i = 0; i < kInvCount; ++i) {
      const Cell& c = table[m.name][i];
      char buf[32];
      if (c.checked == 0) {
        std::snprintf(buf, sizeof buf, "%s", "-");
      } else if (expected(m, static_cast<Inv>(i))) {
        std::snprintf(buf, sizeof buf, "%ld/%ld", c.violated, c.checked);
        realViolations += c.violated;
      } else {
        std::snprintf(buf, sizeof buf, "(%ld/%ld)", c.violated, c.checked);
      }
      std::printf(" %14s", buf);
    }
    std::printf("\n");
  }

  std::printf("\n  violations of invariants the method actually claims: %ld\n", realViolations);
  return realViolations == 0 ? 0 : 1;
}
