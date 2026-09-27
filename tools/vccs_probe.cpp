// VCCS seeding probe: how many seeds each track places, and whether the
// optimized track's seed pruning actually removes anything.
//
// Task 15 measured vccs_optimized over-segmenting E1 by 4.956-6.296x where the
// fair track gives 1.088-1.341x. The candidate explanation on record is that
// adaptive seeding seeds every cell containing structure while the fair track
// seeds only cells whose CENTRE VOXEL happens to be occupied. If that is the
// whole story the cluster count should follow from the seed count, which is
// what this measures.
//
// It runs the shipped algorithms rather than reimplementing their seeding, so
// the numbers cannot drift from the code that produced the campaign. The
// counters it prints are filled in cluster::vccs and cluster::vccsOptimized
// and are inert everywhere else.
//
//   bls_vccs_probe <system.pdb> <config.in> [seed_resolution]
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "config/Parser.hpp"
#include "grid/Grid.hpp"
#include "grid/GridSpec.hpp"
#include "io/TrajectoryReader.hpp"

using namespace bls;

int main(int argc, char** argv) {
  if (argc < 3 || argc > 4) {
    std::fprintf(stderr, "usage: %s <system.pdb> <config.in> [seed_resolution]\n",
                 argv[0]);
    return 2;
  }
  const std::string sys = argv[1], conf = argv[2];
  const double seedRes = (argc == 4) ? std::atof(argv[3]) : 3.0;

  BLSConfig config;
  std::string err;
  Parser parser;
  if (!parser.parseFile(conf, config, err)) {
    std::fprintf(stderr, "config: %s\n", err.c_str());
    return 1;
  }

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

  // The grid of the shared builder, exactly as bls_analyze builds it for every method
  // (audit 27.09.26, finding HB-1: this tool used to fit its own BOX AUTO box, ignoring
  // BOX MANUAL / BOX CELL and the RAM guard; under BOX AUTO the two rules agree).
  GridSpec spec;
  if (!deriveGridSpec(config, frame, nullptr, spec, err)) {
    std::fprintf(stderr, "grid: %s\n", err.c_str());
    return 1;
  }
  const int nx = spec.nx, ny = spec.ny, nz = spec.nz;
  const bool periodic = spec.pbc.all();
  Grid grid;
  configureGrid(grid, spec);
  grid.rasterize(frame.xyz, nullptr, config.cutoff, config.occupancy);

  std::size_t occ = 0;
  for (auto v : grid.occupancy()) occ += (v == 1);

  ClusterParams params;
  params.nx = nx; params.ny = ny; params.nz = nz;
  params.eps = seedRes;
  params.connectivity = config.connectivity;
  params.periodic = periodic;  // deck PBC xyz: same periodicity as the grid and as BLS

  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited());
  // The textbook (fair) track has no periodic form and refuses PBC xyz; it is withdrawn
  // (D16) and reported only under PBC none.
  ClusterResult fair;
  if (!periodic)
    fair = runClusterAlgorithm(ClusterAlgorithm::VCCS, params, grid.occupancy(), grid.visited());
  ClusterResult opt = runClusterAlgorithm(ClusterAlgorithm::VCCSOptimized, params,
                                          grid.occupancy(), grid.visited());

  std::printf("grid %dx%dx%d  volume %zu  occupied %zu  seed_resolution %.3g\n",
              nx, ny, nz, (std::size_t)nx * ny * nz, occ, seedRes);
  std::printf("dfs        nclusters %7d  max_cluster %7d  visited %8zu\n",
              dfs.nclusters, dfs.maxCluster, dfs.visitedVoxels);
  if (periodic)
    std::printf("fair       n/a (textbook track: no periodic form, withdrawn)\n");
  else
    std::printf("fair       nclusters %7d  max_cluster %7d  visited %8zu  "
                "candidates %8d  seeds %7d  prune_min %.4g\n",
                fair.nclusters, fair.maxCluster, fair.visitedVoxels,
                fair.seedCandidates, fair.seedsPlaced, fair.seedPruneThreshold);
  std::printf("optimized  nclusters %7d  max_cluster %7d  visited %8zu  "
              "candidates %8d  seeds %7d  prune_min %.4g\n",
              opt.nclusters, opt.maxCluster, opt.visitedVoxels,
              opt.seedCandidates, opt.seedsPlaced, opt.seedPruneThreshold);
  std::printf("ratios     fair/dfs %.4f  optimized/dfs %.4f  "
              "fair_clusters/seeds %.4f  optimized_clusters/seeds %.4f\n",
              dfs.nclusters ? (double)fair.nclusters / dfs.nclusters : 0.0,
              dfs.nclusters ? (double)opt.nclusters / dfs.nclusters : 0.0,
              fair.seedsPlaced ? (double)fair.nclusters / fair.seedsPlaced : 0.0,
              opt.seedsPlaced ? (double)opt.nclusters / opt.seedsPlaced : 0.0);
  return 0;
}
