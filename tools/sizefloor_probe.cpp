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
// ORIGIN OFFSETS (audit 27.09.26, finding S2-3). An offset is a translation of the
// probe lattice, deck LATTICE_ORIGIN, never of the atoms: the occupancy grid and the
// DFS ground truth are then identical at every offset by construction, and only the
// lattice phase moves. The offset set is named on the command line (tools/offset_set.hpp):
// halton64 = (0,0,0) + 63 Halton points spread over one conventional cell of the probe
// lattice (3-D, sub-voxel); diagonal8 = the pre-audit (k,k,k) set, k = 0..7, expressed
// as the equivalent lattice translation (-k,-k,-k), for the PBC-none reproduction.
//
// THE CELL. PBC none (the pre-audit mechanism, Task 10B): a FIXED cell, supplied the way
// a CRYST1 record supplies one, with the atoms translated inside it. Under BOX AUTO a
// translation of the coordinates alone is a no-op (the origin is minPos), which is why a
// fixed cell is used; its size, including the room the old atom offsets needed, is kept
// so the grid is the pre-audit one. PBC xyz: the deck's cell (BOX CELL / BOX MANUAL),
// untouched. In both cases the DFS grid is built with the shared builder
// (deriveGridSpec, configureGrid) from the same frame BLS gets.
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
//   bls_sizefloor <system.pdb> <config.in> <lattice> <P|I|F> <offset>
//                 --offsets <halton64|diagonal8> [--dump-pdb <out.pdb>]
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
#include "grid/GridSpec.hpp"
#include "io/TrajectoryReader.hpp"
#include "offset_set.hpp"

using namespace bls;

// PBC none: voxels of clearance kept between the atoms and each face of the fixed cell,
// and the room the pre-audit atom offsets (up to 7 voxels) needed. Both stay in the cell
// size so the grid is the pre-audit one; the atoms are no longer moved.
static const int kPadVoxels = 10;
static const int kMaxOffset = 7;

int main(int argc, char** argv) {
  const char* usage =
      "usage: %s <system.pdb> <config.in> <cubic|hexagonal|triclinic> <P|I|F> <offset> "
      "--offsets <halton64|diagonal8> [--dump-pdb <out.pdb>]\n";
  if (argc < 6) {
    std::fprintf(stderr, usage, argv[0]);
    return 2;
  }
  const std::string sys = argv[1], conf = argv[2];
  const std::string latName = argv[3], cenName = argv[4];
  const int offset = std::atoi(argv[5]);
  std::string dumpPdb, offsetSet;
  for (int i = 6; i < argc; ++i) {
    const std::string opt = argv[i];
    if (i + 1 >= argc) { std::fprintf(stderr, usage, argv[0]); return 2; }
    if (opt == "--dump-pdb") dumpPdb = argv[++i];
    else if (opt == "--offsets") offsetSet = argv[++i];
    else { std::fprintf(stderr, "unknown option '%s'\n", opt.c_str()); return 2; }
  }
  const int nOffsets = bls_tools::offsetCount(offsetSet);
  if (nOffsets == 0) {
    std::fprintf(stderr, "--offsets must be halton64 or diagonal8\n");
    return 2;
  }
  if (offset < 0 || offset >= nOffsets) {
    std::fprintf(stderr, "offset must be 0..%d for %s\n", nOffsets - 1, offsetSet.c_str());
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

  const bool periodic = config.pbc.all();
  if (!periodic) {
    // ---- PBC none: the fixed cell, and the atoms translated inside it -------
    const double gs = config.gridSpacing;
    const double pad = kPadVoxels * gs;
    Vec3 mn = frame.xyz[0], mx = frame.xyz[0];
    for (const auto& p : frame.xyz) {
      mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
      mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
    }
    // Rounded UP to a whole number of grid spacings on every axis. Grid sizes its
    // voxels as box/n with n = ceil(box/spacing), so unless the cell is an exact
    // multiple of the spacing the voxel is slightly SMALLER than GRID_SPACING. The
    // pre-audit atom offsets needed the voxel to be exactly GRID_SPACING (measured
    // before this rounding: the occupied-voxel count drifted across offsets, 21416,
    // 21416, 21408, 21399, ...); it is kept so the grid is unchanged.
    auto upToSpacing = [gs](double v) { return std::ceil(v / gs) * gs; };
    const Vec3 cell{upToSpacing(mx.x - mn.x + 2 * pad + kMaxOffset * gs),
                    upToSpacing(mx.y - mn.y + 2 * pad + kMaxOffset * gs),
                    upToSpacing(mx.z - mn.z + 2 * pad + kMaxOffset * gs)};
    const Vec3 shift{pad - mn.x, pad - mn.y, pad - mn.z};
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
  } else if (!dumpPdb.empty()) {
    std::fprintf(stderr, "--dump-pdb is for PBC none (the fixed-cell round trip)\n");
    return 2;
  }

  // ---- the grid, from the shared builder; the offset in the probe basis ----
  GridSpec spec;
  if (!deriveGridSpec(config, frame, nullptr, spec, err)) {
    std::fprintf(stderr, "grid: %s\n", err.c_str());
    return 1;
  }
  const int nx = spec.nx, ny = spec.ny, nz = spec.nz;
  Mat3 basis;
  double dnnVoxel = 0.0;
  if (!Analyzer(config).probeBasis(nx, ny, nz, periodic, basis, dnnVoxel, err)) {
    std::fprintf(stderr, "bls: %s\n", err.c_str());
    return 1;
  }
  const Vec3 origin = bls_tools::latticeOrigin(offsetSet, offset, basis);
  config.latticeOrigin[0] = origin.x;
  config.latticeOrigin[1] = origin.y;
  config.latticeOrigin[2] = origin.z;

  // ---- BLS, through the shipped Analyzer, with labels ---------------------
  Analyzer analyzer(config);
  FrameMetrics metrics;
  std::vector<int> blsLabels;
  if (!analyzer.processFrame(frame, metrics, err, &blsLabels)) {
    std::fprintf(stderr, "bls: %s\n", err.c_str());
    return 1;
  }

  // ---- DFS ground truth, on the same grid spec ----------------------------
  Grid grid;
  configureGrid(grid, spec);
  grid.rasterize(frame.xyz, nullptr, config.cutoff, config.occupancy);

  ClusterParams params;
  params.nx = nx; params.ny = ny; params.nz = nz;
  params.connectivity = config.connectivity;
  params.periodic = periodic;  // deck PBC xyz: same periodicity as the grid and as BLS
  std::vector<int> dfsLabels;
  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited(), &dfsLabels);

  // ---- grid agreement between the BLS path and the DFS path --------------
  // Both paths use the shared builder on the same frame, so the agreement should
  // hold by construction; it is still checked, not assumed.
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

  std::printf("CELL system=%s lattice=%s centering=%s offset=%d offset_set=%s "
              "origin=%.6f,%.6f,%.6f pbc=%s dnn_vox=%.6f "
              "nx=%d ny=%d nz=%d occupied=%zu "
              "dfs_nclusters=%d dfs_max=%d "
              "bls_nclusters=%d bls_max=%d bls_seeds=%d bls_seedhits=%d "
              "present=%d found=%d missed=%d largest_missed=%d smallest_found=%d "
              "truncated=%d split=%d merged=%d grid_ok=%d\n",
              sys.c_str(), latName.c_str(), cenName.c_str(), offset, offsetSet.c_str(),
              origin.x, origin.y, origin.z, periodic ? "xyz" : "none", metrics.dnnVoxel,
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
