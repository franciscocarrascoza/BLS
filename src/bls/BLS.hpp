#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "io/TrajectoryReader.hpp"

#include <memory>

namespace bls {

class Grid;

struct FrameMetrics {
  std::size_t frameIndex{0};
  double timePs{0.0};
  int natoms{0};
  int nx{0}, ny{0}, nz{0};
  double dnnVoxel{0.0};
  std::string lattice;
  std::string centering;
  int seeds{0};
  int seedHits{0};
  int nclusters{0};
  int maxCluster{0};
  std::size_t refinedVoxels{0};
  // Wall-clock of the labelling call alone (manuscript §3.2): for BLS the probe
  // evaluation + refinement, for a comparison method its runClusterAlgorithm call. The
  // grid (box, allocation, voxelisation) is built before this timer starts, identically
  // for every method. Until the 2026-09-26 audit this column timed the whole frame.
  double elapsedMs{0.0};
  // Whole-frame wall-clock: box derivation + grid configure + voxelisation + labelling
  // (the pre-audit scope of elapsed_ms). Fingerprint hashing is outside both timers.
  double totalMs{0.0};
  // BLS only: lattice probe sites evaluated (m of §2.1.2). 0 for other methods.
  long long probes{0};
  // BLS only (audit 28.09.26): elapsedMs split by stage, from one clock read at each boundary.
  //   clearMs  -- zeroing the visited array (every method does this inside its timer);
  //   probeMs  -- stage 1, probe evaluation: site -> voxel -> occupancy, seed list sorted and
  //               de-duplicated (the Enumerator constructor);
  //   refineMs -- stage 2, the refinement walks from the seeds.
  // elapsedMs - (clearMs + probeMs + refineMs) is the basis set-up and the final size sort.
  // Negative (unset) for the comparison methods, which have no such stages.
  double clearMs{-1.0}, probeMs{-1.0}, refineMs{-1.0};
  std::string refinement;  // BLS only: "skip_dfs" or "dfs" (deck REFINEMENT); empty otherwise
  // Grid fingerprint (brief §5): per-axis voxel edge, origin, periodic flags, SHA-256 of
  // the occupancy bits. Identical across methods for the same input and deck.
  double hx{0.0}, hy{0.0}, hz{0.0};
  double originX{0.0}, originY{0.0}, originZ{0.0};
  std::string pbc{"none"};
  std::string occSha256;
  std::vector<int> clusterSizes;
};

class Analyzer {
 public:
  explicit Analyzer(const BLSConfig& config);
  ~Analyzer();

  void setSelection(const std::vector<int>& indices, int natoms);

  // When `labels` is non-null it is resized to nx*ny*nz and overwritten:
  // unoccupied voxels get -1, occupied voxels get dense ids 0..nclusters-1.
  // This is BLS's side of the label contract documented in cluster/Algorithms.hpp.
  // Null (the default) allocates nothing and leaves the timed path untouched.
  // Builds the grid with the shared builder (deriveGridSpec + voxelisation) and labels
  // it; elapsedMs = labelling only, totalMs = whole frame. Used by tools and tests;
  // bls_analyze builds the grid itself and calls labelGrid() directly, exactly as it
  // does for every comparison method.
  bool processFrame(const Frame& frame, FrameMetrics& metrics, std::string& err,
                    std::vector<int>* labels = nullptr);

  // BLS proper on an already voxelised grid: probe evaluation + refinement, timed as
  // metrics.elapsedMs. Initialises the grid's visited array itself, inside the timer,
  // as every comparison method does. Fills seeds, seedHits, nclusters, maxCluster,
  // refinedVoxels, probes, clusterSizes, dnnVoxel, lattice, centering, elapsedMs.
  bool labelGrid(Grid& grid, FrameMetrics& metrics, std::string& err,
                 std::vector<int>* labels = nullptr);

  // The probe-lattice basis (columns, voxels) and dNN (voxels) labelGrid uses on an
  // nx x ny x nz grid: the scaled basis, or under PBC its commensurate cubic form (audit
  // D10; false with err set where labelGrid would refuse the lattice). labelGrid calls
  // this; tools use it to place LATTICE_ORIGIN offsets inside one conventional cell.
  bool probeBasis(int nx, int ny, int nz, bool periodic, Mat3& basis, double& dnnVoxel,
                  std::string& err) const;

  double gridSpacing() const { return config_.gridSpacing; }

 private:
  BLSConfig config_;
  std::vector<int> selection_;
  bool selectionIsAll_{true};

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace bls
