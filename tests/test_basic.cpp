#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "bls/BLS.hpp"
#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "grid/Grid.hpp"
#include "lattice/Basis.hpp"
#include "refine/SkipDFS.hpp"

using bls::Analyzer;
using bls::BLSConfig;
using bls::CenteringType;
using bls::Frame;
using bls::FrameMetrics;
using bls::Grid;
using bls::LatticeDescriptor;
using bls::LatticeSettings;
using bls::LatticeType;
using bls::Mat3;
using bls::SkipDFS;
using bls::SkipDFSConfig;
using bls::Vec3;

namespace {

double computeMinDistance(const Mat3& basis, const std::vector<Vec3>& offsets) {
  double dmin = 1e9;
  for (std::size_t i = 0; i < offsets.size(); ++i) {
    for (int ix = -1; ix <= 1; ++ix) {
      for (int iy = -1; iy <= 1; ++iy) {
        for (int iz = -1; iz <= 1; ++iz) {
          for (std::size_t j = 0; j < offsets.size(); ++j) {
            if (ix == 0 && iy == 0 && iz == 0 && i == j) continue;
            Vec3 diffIdx = offsets[j] + Vec3{static_cast<double>(ix), static_cast<double>(iy),
                                             static_cast<double>(iz)} -
                           offsets[i];
            Vec3 diff = basis * diffIdx;
            double dist = std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
            if (dist > 1e-6) dmin = std::min(dmin, dist);
          }
        }
      }
    }
  }
  return dmin;
}

std::size_t gridIndex(int x, int y, int z, int ny, int nz) {
  return static_cast<std::size_t>(x) * ny * nz + static_cast<std::size_t>(y) * nz +
         static_cast<std::size_t>(z);
}

}  // namespace

int main() {
  int failures = 0;
  auto check = [&](bool cond, const std::string& msg) {
    if (!cond) {
      std::cerr << "FAIL: " << msg << "\n";
      ++failures;
    }
  };

  // Geometry check: FCC lattice scaling.
  {
    LatticeSettings settings;
    settings.lattice = LatticeType::Cubic;
    settings.centering = CenteringType::F;
    LatticeDescriptor desc = bls::buildLattice(settings);
    double scaledTarget = 5.0;
    Mat3 scaled = desc.basis * (scaledTarget / desc.dmin);
    double dminScaled = computeMinDistance(scaled, desc.offsets);
    check(std::abs(dminScaled - scaledTarget) < 1e-3,
          "FCC nearest neighbour mismatch: " + std::to_string(dminScaled));
  }

  // Connectivity test.
  {
    Grid grid;
    grid.configure(2, 2, 1, 1.0, Mat3{Vec3{2, 0, 0}, Vec3{0, 2, 0}, Vec3{0, 0, 1}},
                   Vec3{0, 0, 0});
    auto& occ = grid.occupancy();
    auto& visited = grid.visited();
    occ.assign(grid.size(), 0);
    visited.assign(grid.size(), 0);
    std::size_t idxA = gridIndex(0, 0, 0, 2, 1);
    std::size_t idxB = gridIndex(1, 1, 0, 2, 1);
    occ[idxA] = 1;
    occ[idxB] = 1;

    SkipDFSConfig cfg6{2, 2, 1, 6, 1};
    SkipDFS dfs6(cfg6, occ, visited);
    int size6 = dfs6.runFrom(0, 0, 0);
    check(size6 == 1, "Connectivity 6 should not connect diagonal voxels");

    visited.assign(grid.size(), 0);
    SkipDFSConfig cfg26{2, 2, 1, 26, 1};
    SkipDFS dfs26(cfg26, occ, visited);
    int size26 = dfs26.runFrom(0, 0, 0);
    check(size26 == 2, "Connectivity 26 should connect diagonal voxels");
  }

  // PBC wrapping test.
  //
  // This test used to place atoms at x = -0.1 AND x = 3.9 and assert only that
  // voxel 3 ended up occupied. Both atoms land in voxel 3 -- the second one
  // without any wrapping at all -- so the assertion held whether or not
  // wrapping worked, and it kept holding after wrapping was made conditional
  // on BoxPeriodicity. The out-of-box atom is now tested on its own, and the
  // non-periodic case is tested for the opposite outcome.
  {
    const Mat3 box{Vec3{4, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    const std::vector<Vec3> outside = {Vec3{-0.1, 0.0, 0.0}};

    Grid periodic;
    periodic.configure(4, 1, 1, 1.0, box, Vec3{0, 0, 0}, bls::BoxPeriodicity::Periodic);
    periodic.rasterize(outside, nullptr, 0.0, bls::OccupancyMode::Any);
    check(periodic.occupancy()[gridIndex(3, 0, 0, 1, 1)] == 1,
          "PBC wrapping failed to map an out-of-box atom into the far voxel");

    Grid nonPeriodic;
    nonPeriodic.configure(4, 1, 1, 1.0, box, Vec3{0, 0, 0}, bls::BoxPeriodicity::NonPeriodic);
    nonPeriodic.rasterize(outside, nullptr, 0.0, bls::OccupancyMode::Any);
    std::size_t occupiedCount = 0;
    for (uint8_t v : nonPeriodic.occupancy())
      if (v) ++occupiedCount;
    check(occupiedCount == 0,
          "A non-periodic box must drop an out-of-box atom, not fold it onto a face");

    // An in-box atom is unaffected by either setting.
    for (Grid* g : {&periodic, &nonPeriodic}) {
      g->rasterize({Vec3{3.9, 0.0, 0.0}}, nullptr, 0.0, bls::OccupancyMode::Any);
      check(g->occupancy()[gridIndex(3, 0, 0, 1, 1)] == 1,
            "In-box atom was not rasterised into its own voxel");
    }
  }

  // Determinism test.
  {
    BLSConfig cfg;
    cfg.gridSpacing = 1.0;
    cfg.connectivity = 6;
    cfg.refinementStride = 1;
    cfg.lattice.lattice = LatticeType::Cubic;
    cfg.lattice.centering = CenteringType::P;

    Analyzer analyzer(cfg);

    Frame frame;
    frame.xyz = {Vec3{0.1, 0.0, 0.0}, Vec3{1.1, 0.0, 0.0}};
    frame.natoms = static_cast<int>(frame.xyz.size());
    frame.time = 0.0;
    frame.box = Mat3{Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}, Vec3{0.0, 0.0, 2.0}};

    analyzer.setSelection({}, frame.natoms);

    FrameMetrics m1;
    FrameMetrics m2;
    std::string err;
    check(analyzer.processFrame(frame, m1, err), "First deterministic run failed: " + err);
    check(analyzer.processFrame(frame, m2, err), "Second deterministic run failed: " + err);
    check(m1.nclusters == m2.nclusters && m1.maxCluster == m2.maxCluster,
          "Determinism violated: metrics differ");
  }

  // CC3D connectivity test.
  {
    std::vector<uint8_t> occupancy(8, 0);
    std::vector<uint8_t> visited(8, 0);
    // Create two voxels at opposite corners
    occupancy[gridIndex(0, 0, 0, 2, 2)] = 1;
    occupancy[gridIndex(1, 1, 1, 2, 2)] = 1;

    // Test 6-connectivity (face neighbors only)
    bls::ClusterResult result6 = bls::cc3d(2, 2, 2, 6, occupancy, visited);
    check(result6.nclusters == 2, "CC3D with 6-connectivity should find 2 separate clusters");

    // Reset visited
    visited.assign(8, 0);

    // Test 26-connectivity (includes diagonal)
    bls::ClusterResult result26 = bls::cc3d(2, 2, 2, 26, occupancy, visited);
    check(result26.nclusters == 1, "CC3D with 26-connectivity should find 1 cluster");
  }

  // HDBSCAN basic functionality test.
  {
    std::vector<uint8_t> occupancy(27, 0);
    std::vector<uint8_t> visited(27, 0);
    // Create a small cluster in the center
    occupancy[gridIndex(1, 1, 1, 3, 3)] = 1;
    occupancy[gridIndex(1, 1, 0, 3, 3)] = 1;
    occupancy[gridIndex(1, 1, 2, 3, 3)] = 1;
    occupancy[gridIndex(1, 0, 1, 3, 3)] = 1;
    occupancy[gridIndex(1, 2, 1, 3, 3)] = 1;

    bls::ClusterResult result = bls::hdbscan(3, 3, 3, 2, 2, occupancy, visited);
    check(result.nclusters >= 1, "HDBSCAN should find at least 1 cluster");
    check(result.visitedVoxels > 0, "HDBSCAN should visit some voxels");
  }

  // Algorithm comparison test (GCBD vs CC3D should match with same connectivity).
  {
    std::vector<uint8_t> occupancy(8, 0);
    std::vector<uint8_t> visited1(8, 0);
    std::vector<uint8_t> visited2(8, 0);
    // Create a connected pattern
    occupancy[gridIndex(0, 0, 0, 2, 2)] = 1;
    occupancy[gridIndex(0, 0, 1, 2, 2)] = 1;
    occupancy[gridIndex(0, 1, 0, 2, 2)] = 1;

    bls::ClusterResult gcbd = bls::gcbd(2, 2, 2, occupancy, visited1);
    bls::ClusterResult cc3d = bls::cc3d(2, 2, 2, 6, occupancy, visited2);

    check(gcbd.nclusters == cc3d.nclusters,
          "GCBD and CC3D with 6-connectivity should produce same cluster count");
    check(gcbd.maxCluster == cc3d.maxCluster,
          "GCBD and CC3D should produce same max cluster size");
  }

  if (failures == 0) {
    std::cout << "All tests passed.\n";
  }
  return failures == 0 ? 0 : 1;
}
