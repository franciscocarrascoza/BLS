// E7 probe: recall and timing for BLS and traditional_dfs at an arbitrary
// (LATTICE, CENTERING, ALPHA, GRID_SPACING), over the Task 10B/19 origin-
// offset mechanism.
//
// WHY THIS EXISTS. Task E7 needs "recall == 1.0 on all four E1 systems, worst
// case over 8 origin offsets" as a function of ALPHA and GRID_SPACING. No
// existing tool sweeps ALPHA or GRID_SPACING independently of the campaign
// decks, and none reports both BLS and DFS timing from the SAME grid build so
// a ratio is meaningful. This generalises tools/sizefloor_probe.cpp (which
// fixes ALPHA/GRID_SPACING and sweeps LATTICE/CENTERING) to instead accept
// all four as CLI overrides, keeping only CUTOFF, CONNECTIVITY, SKIP and
// OCCUPANCY from the deck.
//
// GRID AGREEMENT AND TIMING CONVENTION. Both algorithms are timed to match
// main.cpp exactly: BLS's Analyzer::processFrame times box setup + rasterize
// + enumeration + refinement (nothing else wraps it); the DFS path here wraps
// a ScopedTimer around Grid::configure + Grid::rasterize + runClusterAlgorithm,
// which is the same span main.cpp's non-BLS branch times (its outer `timer`
// starts before grid setup, per main.cpp:427-566). Without matching that
// convention a ratio_dfs computed from this probe would not be comparable to
// one from the campaign.
//
// ORIGIN OFFSETS: identical mechanism to tools/sizefloor_probe.cpp -- see its
// header for the reasoning (BOX AUTO makes coordinate translation alone a
// no-op; a fixed cell rounded up to a whole number of GRID_SPACING is
// required so a k*GRID_SPACING translation is exactly k voxels).
//
//   bls_e7 <system.pdb> <config.in> <lattice> <P|I|F> <alpha> <grid_spacing>
//          <offset 0..7>
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
#include "io/TrajectoryReader.hpp"
#include "util/Timer.hpp"

using namespace bls;

static const int kPadVoxels = 10;
static const int kMaxOffset = 7;

int main(int argc, char** argv) {
  if (argc != 8) {
    std::fprintf(stderr,
                 "usage: %s <system.pdb> <config.in> "
                 "<cubic|hexagonal|triclinic> <P|I|F> <alpha> "
                 "<grid_spacing> <offset 0..7>\n", argv[0]);
    return 2;
  }
  const std::string sys = argv[1], conf = argv[2];
  const std::string latName = argv[3], cenName = argv[4];
  const double alpha = std::atof(argv[5]);
  const double gs = std::atof(argv[6]);
  const int offset = std::atoi(argv[7]);
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

  // ---- the fixed cell, and the atoms translated inside it -----------------
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
  const Vec3 shift{pad - mn.x + offset * gs,
                   pad - mn.y + offset * gs,
                   pad - mn.z + offset * gs};
  for (auto& p : frame.xyz) p += shift;
  frame.box = Mat3{Vec3{cell.x, 0, 0}, Vec3{0, cell.y, 0}, Vec3{0, 0, cell.z}};

  if (cell.x > (mx.x - mn.x) * 10.0 || cell.y > (mx.y - mn.y) * 10.0 ||
      cell.z > (mx.z - mn.z) * 10.0) {
    std::fprintf(stderr, "fixed cell exceeds 10x the coordinate extent; "
                         "BOX AUTO would discard it\n");
    return 1;
  }

  // ---- BLS, through the shipped Analyzer -----------------------------------
  Analyzer analyzer(config);
  FrameMetrics metrics;
  if (!analyzer.processFrame(frame, metrics, err)) {
    std::fprintf(stderr, "bls: %s\n", err.c_str());
    return 1;
  }

  // ---- DFS ground truth, timed the way main.cpp times the non-BLS path ----
  ScopedTimer dfsTimer;
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
  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited());
  const double dfsElapsedMs = dfsTimer.elapsedMilliseconds();

  int gridOk = (metrics.nx == nx && metrics.ny == ny && metrics.nz == nz) ? 1 : 0;
  std::size_t occupied = 0;
  for (auto v : grid.occupancy()) occupied += (v != 0);

  const int present = dfs.nclusters;
  const int detected = metrics.nclusters;
  const double recall = present > 0 ? (double)detected / present : -1.0;

  std::printf("CELL system=%s lattice=%s centering=%s alpha=%.4g "
              "grid_spacing=%.4g offset=%d nx=%d ny=%d nz=%d occupied=%zu "
              "present=%d detected=%d recall=%.6f "
              "bls_max=%d bls_seeds=%d bls_seedhits=%d bls_elapsed_ms=%.4f "
              "dfs_max=%d dfs_elapsed_ms=%.4f grid_ok=%d\n",
              sys.c_str(), latName.c_str(), cenName.c_str(), alpha, gs, offset,
              nx, ny, nz, occupied, present, detected, recall,
              metrics.maxCluster, metrics.seeds, metrics.seedHits, metrics.elapsedMs,
              dfs.maxCluster, dfsElapsedMs, gridOk);
  return 0;
}
