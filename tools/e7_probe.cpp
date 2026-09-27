// E7 probe: recall and timing for BLS and traditional_dfs at an arbitrary
// (LATTICE, CENTERING, ALPHA, GRID_SPACING), over the Task 10B/19 origin-
// offset mechanism.
//
// WHY THIS EXISTS. Task E7 needs "recall == 1.0 on all four E1 systems, worst
// case over the origin offsets" as a function of ALPHA and GRID_SPACING. No
// existing tool sweeps ALPHA or GRID_SPACING independently of the campaign
// decks, and none reports both BLS and DFS timing from the SAME grid build so
// a ratio is meaningful. This generalises tools/sizefloor_probe.cpp (which
// fixes ALPHA/GRID_SPACING and sweeps LATTICE/CENTERING) to instead accept
// all four as CLI overrides, keeping only CUTOFF, CONNECTIVITY, SKIP and
// OCCUPANCY from the deck.
//
// GRID AND TIMING CONVENTION (audit 27.09.26). Both labellings run on the grid of the
// shared builder (deriveGridSpec), from the same frame. Timing follows manuscript §3.2
// and audit decision D9, as bls_analyze does: bls_elapsed_ms = Analyzer's labelling call
// (probe evaluation + refinement), dfs_elapsed_ms = the runClusterAlgorithm call alone.
// The whole-frame spans (grid derivation + voxelisation + labelling) are printed as
// bls_total_ms / dfs_total_ms. Until this change dfs_elapsed_ms still included the grid
// build while bls_elapsed_ms no longer did (since D9, inner e5ba8c2), so a ratio from this
// probe would not have been comparable (finding HARN-10; no run used it).
//
// ORIGIN OFFSETS AND THE CELL: identical to tools/sizefloor_probe.cpp -- see its header.
// Offsets translate the probe lattice (LATTICE_ORIGIN), the set is named with --offsets
// (tools/offset_set.hpp); PBC none keeps the pre-audit fixed cell, PBC xyz the deck's cell.
//
//   bls_e7 <system.pdb> <config.in> <lattice> <P|I|F> <alpha> <grid_spacing>
//          <offset> --offsets <halton64|diagonal8>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
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
#include "util/Timer.hpp"

using namespace bls;

// PBC none: clearance and the pre-audit offset room, kept in the cell size (see
// tools/sizefloor_probe.cpp).
static const int kPadVoxels = 10;
static const int kMaxOffset = 7;

int main(int argc, char** argv) {
  if (argc != 10 || std::string(argv[8]) != "--offsets") {
    std::fprintf(stderr,
                 "usage: %s <system.pdb> <config.in> "
                 "<cubic|hexagonal|triclinic> <P|I|F> <alpha> "
                 "<grid_spacing> <offset> --offsets <halton64|diagonal8>\n", argv[0]);
    return 2;
  }
  const std::string sys = argv[1], conf = argv[2];
  const std::string latName = argv[3], cenName = argv[4];
  const double alpha = std::atof(argv[5]);
  const double gs = std::atof(argv[6]);
  const int offset = std::atoi(argv[7]);
  const std::string offsetSet = argv[9];
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
  config.alpha = alpha;
  config.gridSpacing = gs;

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
    const double pad = kPadVoxels * gs;
    Vec3 mn = frame.xyz[0], mx = frame.xyz[0];
    for (const auto& p : frame.xyz) {
      mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
      mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
    }
    auto upToSpacing = [gs](double v) { return std::ceil(v / gs) * gs; };
    const Vec3 cell{upToSpacing(mx.x - mn.x + 2 * pad + kMaxOffset * gs),
                    upToSpacing(mx.y - mn.y + 2 * pad + kMaxOffset * gs),
                    upToSpacing(mx.z - mn.z + 2 * pad + kMaxOffset * gs)};
    const Vec3 shift{pad - mn.x, pad - mn.y, pad - mn.z};
    for (auto& p : frame.xyz) p += shift;
    frame.box = Mat3{Vec3{cell.x, 0, 0}, Vec3{0, cell.y, 0}, Vec3{0, 0, cell.z}};

    if (cell.x > (mx.x - mn.x) * 10.0 || cell.y > (mx.y - mn.y) * 10.0 ||
        cell.z > (mx.z - mn.z) * 10.0) {
      std::fprintf(stderr, "fixed cell exceeds 10x the coordinate extent; "
                           "BOX AUTO would discard it\n");
      return 1;
    }
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

  // ---- BLS, through the shipped Analyzer -----------------------------------
  Analyzer analyzer(config);
  FrameMetrics metrics;
  if (!analyzer.processFrame(frame, metrics, err)) {
    std::fprintf(stderr, "bls: %s\n", err.c_str());
    return 1;
  }

  // ---- DFS ground truth: labelling timed alone, whole frame as total -------
  ScopedTimer dfsTotalTimer;
  GridSpec dfsSpec;
  if (!deriveGridSpec(config, frame, nullptr, dfsSpec, err)) {
    std::fprintf(stderr, "grid: %s\n", err.c_str());
    return 1;
  }
  Grid grid;
  configureGrid(grid, dfsSpec);
  grid.rasterize(frame.xyz, nullptr, config.cutoff, config.occupancy);

  ClusterParams params;
  params.nx = nx; params.ny = ny; params.nz = nz;
  params.connectivity = config.connectivity;
  params.periodic = periodic;  // deck PBC xyz: same periodicity as the grid and as BLS
  ScopedTimer dfsTimer;
  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited());
  const double dfsElapsedMs = dfsTimer.elapsedMilliseconds();
  const double dfsTotalMs = dfsTotalTimer.elapsedMilliseconds();

  int gridOk = (metrics.nx == nx && metrics.ny == ny && metrics.nz == nz) ? 1 : 0;
  std::size_t occupied = 0;
  for (auto v : grid.occupancy()) occupied += (v != 0);

  const int present = dfs.nclusters;
  const int detected = metrics.nclusters;
  const double recall = present > 0 ? (double)detected / present : -1.0;

  std::printf("CELL system=%s lattice=%s centering=%s alpha=%.4g "
              "grid_spacing=%.4g offset=%d offset_set=%s origin=%.6f,%.6f,%.6f pbc=%s "
              "dnn_vox=%.6f nx=%d ny=%d nz=%d occupied=%zu "
              "present=%d detected=%d recall=%.6f "
              "bls_max=%d bls_seeds=%d bls_seedhits=%d bls_elapsed_ms=%.4f "
              "dfs_max=%d dfs_elapsed_ms=%.4f bls_total_ms=%.4f dfs_total_ms=%.4f grid_ok=%d\n",
              sys.c_str(), latName.c_str(), cenName.c_str(), alpha, gs, offset,
              offsetSet.c_str(), origin.x, origin.y, origin.z, periodic ? "xyz" : "none",
              metrics.dnnVoxel, nx, ny, nz, occupied, present, detected, recall,
              metrics.maxCluster, metrics.seeds, metrics.seedHits, metrics.elapsedMs,
              dfs.maxCluster, dfsElapsedMs, metrics.totalMs, dfsTotalMs, gridOk);
  return 0;
}
