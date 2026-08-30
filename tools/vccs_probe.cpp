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

  // Same box handling as main.cpp's comparison path: BOX AUTO fits the
  // coordinates, with the same padding, so the grid matches the campaign's.
  Mat3 box = frame.box;
  Vec3 origin{0.0, 0.0, 0.0};
  auto len = [](const Vec3& v) { return norm(v); };
  bool needFit = len(box.column(0)) < 1e-8 || len(box.column(1)) < 1e-8 ||
                 len(box.column(2)) < 1e-8;
  if (!needFit && config.boxMode == BoxMode::Auto && !frame.xyz.empty()) {
    Vec3 mn = frame.xyz[0], mx = frame.xyz[0];
    for (const auto& p : frame.xyz) {
      mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
      mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
    }
    const double ex = mx.x - mn.x, ey = mx.y - mn.y, ez = mx.z - mn.z;
    if (len(box.column(0)) > ex * 10.0 || len(box.column(1)) > ey * 10.0 ||
        len(box.column(2)) > ez * 10.0)
      needFit = true;
  }
  if (needFit && !frame.xyz.empty()) {
    Vec3 mn = frame.xyz[0], mx = frame.xyz[0];
    for (const auto& p : frame.xyz) {
      mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
      mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
    }
    const double pad = config.gridSpacing * 2.0;
    origin = mn;
    box = Mat3{Vec3{std::max(mx.x - mn.x + pad, pad), 0, 0},
               Vec3{0, std::max(mx.y - mn.y + pad, pad), 0},
               Vec3{0, 0, std::max(mx.z - mn.z + pad, pad)}};
  }

  const int nx = std::max(1, (int)std::ceil(norm(box.column(0)) / config.gridSpacing));
  const int ny = std::max(1, (int)std::ceil(norm(box.column(1)) / config.gridSpacing));
  const int nz = std::max(1, (int)std::ceil(norm(box.column(2)) / config.gridSpacing));

  Grid grid;
  grid.configure(nx, ny, nz, config.gridSpacing, box, origin,
                 periodicityForBoxMode(config.boxMode));
  grid.rasterize(frame.xyz, nullptr, config.cutoff, config.occupancy);

  std::size_t occ = 0;
  for (auto v : grid.occupancy()) occ += (v == 1);

  ClusterParams params;
  params.nx = nx; params.ny = ny; params.nz = nz;
  params.eps = seedRes;
  params.connectivity = config.connectivity;

  ClusterResult dfs = runClusterAlgorithm(ClusterAlgorithm::TraditionalDFS, params,
                                          grid.occupancy(), grid.visited());
  ClusterResult fair = runClusterAlgorithm(ClusterAlgorithm::VCCS, params,
                                           grid.occupancy(), grid.visited());
  ClusterResult opt = runClusterAlgorithm(ClusterAlgorithm::VCCSOptimized, params,
                                          grid.occupancy(), grid.visited());

  std::printf("grid %dx%dx%d  volume %zu  occupied %zu  seed_resolution %.3g\n",
              nx, ny, nz, (std::size_t)nx * ny * nz, occ, seedRes);
  std::printf("dfs        nclusters %7d  max_cluster %7d  visited %8zu\n",
              dfs.nclusters, dfs.maxCluster, dfs.visitedVoxels);
  std::printf("fair       nclusters %7d  max_cluster %7d  visited %8zu  "
              "candidates %8d  seeds %7d  prune_min %d\n",
              fair.nclusters, fair.maxCluster, fair.visitedVoxels,
              fair.seedCandidates, fair.seedsPlaced, fair.seedPruneThreshold);
  std::printf("optimized  nclusters %7d  max_cluster %7d  visited %8zu  "
              "candidates %8d  seeds %7d  prune_min %d\n",
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
