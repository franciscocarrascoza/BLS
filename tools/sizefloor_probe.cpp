// Size-floor probe: for every connected component in the DFS ground truth, its
// size in voxels and whether BLS found it.
//
// WHY THIS EXISTS. The paper's central claim is that BLS returns exact
// connected components ABOVE A STATED SIZE FLOOR. Task 10B measured that floor
// but ran diagnostically and wrote nothing to disk; the campaign CSVs carry
// nclusters and max_cluster and never a per-component size, so the claim has
// existed only as prose. Cluster counts cannot support it: a run that finds 984
// of 991 says nothing about WHICH seven were lost or how big they were.
//
// It runs the shipped Analyzer and the shipped cluster::traditionalDFS through
// their label contract (Algorithms.hpp: unoccupied -1, occupied a dense id), so
// the partition measured here is the partition the campaign produced. Nothing
// is reimplemented.
//
// ORIGIN OFFSETS. BLS's lattice is anchored to the grid origin, so the answer
// depends on where the lattice falls relative to the structure. Translating the
// coordinates alone is a NO-OP under BOX AUTO: with no usable cell the origin
// is minPos, so the origin translates with the atoms and the lattice phase is
// unchanged. BOX MANUAL is defect 12 (main.cpp ignores it). The mechanism used
// here is the Task 10B one: a FIXED cell, supplied the way a CRYST1 record
// supplies one, with the atoms translated inside it. The cell is identical for
// every offset and every lattice, so nx, ny, nz and the voxelisation are fixed
// and only the lattice phase moves.
//
// Offsets are INTEGER voxel translations (k,k,k) for k = 0..7. Two reasons:
//   - an integer-voxel translation is a pure index shift of the occupancy, so
//     the DFS ground truth is identical across offsets and the size histogram
//     has a stable denominator. Anything else would change what is being
//     measured at the same time as the lattice phase.
//   - the lattice period is dNN = 2.4 voxels at the E1 deck, which is not
//     commensurate with the voxel grid, so k mod 2.4 for k = 0..7 samples eight
//     distinct phases spread across the period (0, 1, 2, 0.6, 1.6, 0.2, 1.2,
//     2.2). A {0,1}^3 corner set would sample only two.
//
// Recall and component sizes are integers produced by a deterministic
// traversal of a deterministic grid. They do not vary between runs, so this
// probe takes no replicates; timing is not measured here and is not the point.
//
// EQUIVALENCE. The fixed cell is applied in memory rather than by writing a
// PDB, which avoids the %8.3f coordinate round-trip. --dump-pdb writes exactly
// the cell and coordinates this probe uses, as a CRYST1 record and ATOM lines,
// so the mechanism can be replayed through the shipped bls_analyze and the two
// paths compared. run_sizefloor.sh does that once per sweep.
//
//   bls_sizefloor <system.pdb> <config.in> <lattice> <P|I|F> <offset 0..7>
//                 [--dump-pdb <out.pdb>]
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "bls/BLS.hpp"
#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "config/Parser.hpp"
#include "grid/Grid.hpp"
#include "io/TrajectoryReader.hpp"

using namespace bls;

// Voxels of clearance kept between the atoms and each face of the fixed cell.
// It has to exceed the largest offset (7 voxels) so that no offset pushes an
// atom out of the cell: the analysis box is NonPeriodic under BOX AUTO, so an
// atom outside it would be clipped rather than wrapped, and the ground truth
// would then differ between offsets for a reason that has nothing to do with
// the lattice.
static const int kPadVoxels = 10;
static const int kMaxOffset = 7;

int main(int argc, char** argv) {
  if (argc != 6 && argc != 8) {
    std::fprintf(stderr,
                 "usage: %s <system.pdb> <config.in> "
                 "<cubic|hexagonal|triclinic> <P|I|F> <offset 0..7> "
                 "[--dump-pdb <out.pdb>]\n", argv[0]);
    return 2;
  }
  const std::string sys = argv[1], conf = argv[2];
  const std::string latName = argv[3], cenName = argv[4];
  const int offset = std::atoi(argv[5]);
  std::string dumpPdb;
  if (argc == 8) {
    if (std::string(argv[6]) != "--dump-pdb") {
      std::fprintf(stderr, "unknown option '%s'\n", argv[6]);
      return 2;
    }
    dumpPdb = argv[7];
  }
  if (offset < 0 || offset > kMaxOffset) {
    std::fprintf(stderr, "offset must be 0..%d\n", kMaxOffset);
    return 2;
  }

  BLSConfig config;
  std::string err;
  Parser parser;
  if (!parser.parseFile(conf, config, err)) {
    std::fprintf(stderr, "config: %s\n", err.c_str());
    return 1;
  }
  // Override only the lattice and centering. Everything else -- grid spacing,
  // alpha, cutoff, connectivity, refinement stride, occupancy mode -- stays at
  // the campaign deck's value, so the six combinations differ in nothing else.
  if      (latName == "cubic")      config.lattice.lattice = LatticeType::Cubic;
  else if (latName == "hexagonal")  config.lattice.lattice = LatticeType::Hexagonal;
  else if (latName == "triclinic")  config.lattice.lattice = LatticeType::Triclinic;
  else { std::fprintf(stderr, "unknown lattice '%s'\n", latName.c_str()); return 2; }
  if      (cenName == "P") config.lattice.centering = CenteringType::P;
  else if (cenName == "I") config.lattice.centering = CenteringType::I;
  else if (cenName == "F") config.lattice.centering = CenteringType::F;
  else { std::fprintf(stderr, "unknown centering '%s'\n", cenName.c_str()); return 2; }
  config.lattice.latticeSet = true;
  config.lattice.centeringSet = true;

  auto reader = makeTrajectoryReader(sys, "auto", err);
  if (!reader || !reader->open(sys, err)) {
    std::fprintf(stderr, "system: %s\n", err.c_str());
    return 1;
  }
  Frame frame;
  if (!reader->read(frame, err)) {
    std::fprintf(stderr, "read: %s\n", err.c_str());
    return 1;
  }
  if (frame.xyz.empty()) { std::fprintf(stderr, "empty frame\n"); return 1; }

  // ---- the fixed cell, and the atoms translated inside it -----------------
  const double gs = config.gridSpacing;
  const double pad = kPadVoxels * gs;
  Vec3 mn = frame.xyz[0], mx = frame.xyz[0];
  for (const auto& p : frame.xyz) {
    mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
    mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
  }
  // Cell dimensions do not depend on the offset: extent, clearance on both
  // faces, and room for the largest offset. Every cell of the matrix therefore
  // rasterises onto exactly the same grid.
  // Rounded UP to a whole number of grid spacings on every axis. Grid sizes
  // its voxels as box/n with n = ceil(box/spacing), so unless the cell is an
  // exact multiple of the spacing the voxel is slightly SMALLER than
  // GRID_SPACING and a translation by k*GRID_SPACING is not a translation by k
  // voxels. Measured before this rounding was added: the occupied-voxel count
  // drifted across offsets (21416, 21416, 21408, 21399, ...), which would have
  // moved the ground truth and the lattice phase at the same time. With the
  // cell an exact multiple, box/n == spacing and the shift is exactly k voxels.
  auto upToSpacing = [gs](double v) { return std::ceil(v / gs) * gs; };
  const Vec3 cell{upToSpacing(mx.x - mn.x + 2 * pad + kMaxOffset * gs),
                  upToSpacing(mx.y - mn.y + 2 * pad + kMaxOffset * gs),
                  upToSpacing(mx.z - mn.z + 2 * pad + kMaxOffset * gs)};
  const Vec3 shift{pad - mn.x + offset * gs,
                   pad - mn.y + offset * gs,
                   pad - mn.z + offset * gs};
  for (auto& p : frame.xyz) p += shift;
  frame.box = Mat3{Vec3{cell.x, 0, 0}, Vec3{0, cell.y, 0}, Vec3{0, 0, cell.z}};

  // A cell more than 10x the coordinate extent is rejected by BOX AUTO as
  // implausible and replaced by a fitted bounding box, which would silently
  // undo the whole mechanism. Check rather than assume.
  if (cell.x > (mx.x - mn.x) * 10.0 || cell.y > (mx.y - mn.y) * 10.0 ||
      cell.z > (mx.z - mn.z) * 10.0) {
    std::fprintf(stderr, "fixed cell exceeds 10x the coordinate extent; "
                         "BOX AUTO would discard it\n");
    return 1;
  }

  if (!dumpPdb.empty()) {
    // Written after the translation, so the file carries the same cell and the
    // same coordinates the measurement below uses. Water only: the E1 systems
    // are O and H, and the element column is not read back by PdbReader.
    std::FILE* fh = std::fopen(dumpPdb.c_str(), "w");
    if (!fh) { std::fprintf(stderr, "cannot write %s\n", dumpPdb.c_str()); return 1; }
    std::fprintf(fh, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
                 cell.x, cell.y, cell.z, 90.0, 90.0, 90.0);
    for (std::size_t i = 0; i < frame.xyz.size(); ++i) {
      const Vec3& p = frame.xyz[i];
      std::fprintf(fh,
                   "ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           O\n",
                   (int)((i % 99999) + 1), (int)((i / 3 % 9999) + 1), p.x, p.y, p.z);
    }
    std::fprintf(fh, "END\n");
    std::fclose(fh);
  }

  // ---- BLS, through the shipped Analyzer, with labels ---------------------
  Analyzer analyzer(config);
  FrameMetrics metrics;
  std::vector<int> blsLabels;
  if (!analyzer.processFrame(frame, metrics, err, &blsLabels)) {
    std::fprintf(stderr, "bls: %s\n", err.c_str());
    return 1;
  }

  // ---- DFS ground truth, on a grid built the same way ---------------------
  const int nx = std::max(1, (int)std::ceil(cell.x / gs));
  const int ny = std::max(1, (int)std::ceil(cell.y / gs));
  const int nz = std::max(1, (int)std::ceil(cell.z / gs));
  Grid grid;
  grid.configure(nx, ny, nz, gs, frame.box, Vec3{0, 0, 0},
                 periodicityForBoxMode(config.boxMode));
  grid.rasterize(frame.xyz, nullptr, config.cutoff, config.occupancy);

  ClusterParams params;
  params.nx = nx; params.ny = ny; params.nz = nz;
  params.connectivity = config.connectivity;
  std::vector<int> dfsLabels;
  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited(), &dfsLabels);

  // ---- grid agreement between the BLS path and the DFS path --------------
  // The two paths build their grids from the same box by the same rule but in
  // different translation units, so the agreement is checked, not assumed.
  // Dimensions first, then the stronger test: every voxel BLS labelled must be
  // occupied on the DFS grid. A disagreement in origin, spacing or stencil
  // would put at least one BLS label on a voxel this grid calls empty.
  int gridOk = 1;
  if (metrics.nx != nx || metrics.ny != ny || metrics.nz != nz) gridOk = 0;
  const std::size_t n = (std::size_t)nx * ny * nz;
  std::size_t occupied = 0;
  for (auto v : grid.occupancy()) occupied += (v != 0);
  if (gridOk) {
    for (std::size_t i = 0; i < n; ++i) {
      if (blsLabels[i] >= 0 && dfsLabels[i] < 0) { gridOk = 0; break; }
    }
  }

  // ---- per-component accounting ------------------------------------------
  const int nd = dfs.nclusters;
  std::vector<int> size(nd, 0);          // voxels in this DFS component
  std::vector<int> labelled(nd, 0);      // of those, how many BLS labelled
  std::vector<std::set<int>> blsIds(nd); // which BLS ids appear on it
  std::map<int, std::set<int>> dfsOfBls; // BLS id -> DFS components it spans
  for (std::size_t i = 0; i < n; ++i) {
    const int d = dfsLabels[i];
    if (d < 0) continue;
    ++size[d];
    const int b = blsLabels[i];
    if (b >= 0) {
      ++labelled[d];
      blsIds[d].insert(b);
      dfsOfBls[b].insert(d);
    }
  }

  // A component is FOUND when BLS labelled any of its voxels: BLS reaches a
  // component only by seeding into it, and the size floor is about whether a
  // lattice site lands in the component at all.
  //
  // The three refinement-loss modes are separated because they fail
  // differently. TRUNCATED: BLS entered the component and labelled only part
  // of it -- refinement stopped early. SPLIT: BLS gave one component more than
  // one id -- refinement broke it in two. MERGED: one BLS id covers more than
  // one DFS component -- refinement crossed a gap it should not have.
  std::map<int, long long> presentBySize, missedBySize;
  int found = 0, missed = 0, truncated = 0, split = 0;
  int largestMissed = 0, smallestFound = -1;
  for (int c = 0; c < nd; ++c) {
    if (labelled[c] > 0) {
      ++found;
      if (smallestFound < 0 || size[c] < smallestFound) smallestFound = size[c];
      if (labelled[c] < size[c]) ++truncated;
      if (blsIds[c].size() > 1) ++split;
    } else {
      ++missed;
      largestMissed = std::max(largestMissed, size[c]);
      missedBySize[size[c]]++;
    }
    presentBySize[size[c]]++;
  }
  int merged = 0;
  for (const auto& kv : dfsOfBls) if (kv.second.size() > 1) ++merged;

  std::printf("CELL system=%s lattice=%s centering=%s offset=%d "
              "nx=%d ny=%d nz=%d occupied=%zu "
              "dfs_nclusters=%d dfs_max=%d "
              "bls_nclusters=%d bls_max=%d bls_seeds=%d bls_seedhits=%d "
              "present=%d found=%d missed=%d largest_missed=%d smallest_found=%d "
              "truncated=%d split=%d merged=%d grid_ok=%d\n",
              sys.c_str(), latName.c_str(), cenName.c_str(), offset,
              nx, ny, nz, occupied,
              dfs.nclusters, dfs.maxCluster,
              metrics.nclusters, metrics.maxCluster, metrics.seeds, metrics.seedHits,
              nd, found, missed, largestMissed, smallestFound,
              truncated, split, merged, gridOk);
  for (const auto& kv : presentBySize) {
    long long m = missedBySize.count(kv.first) ? missedBySize[kv.first] : 0;
    std::printf("HIST %d %lld %lld\n", kv.first, kv.second, m);
  }
  return 0;
}
