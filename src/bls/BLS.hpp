#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "io/TrajectoryReader.hpp"

namespace bls {

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
  double elapsedMs{0.0};
  std::vector<int> clusterSizes;
};

class Analyzer {
 public:
  explicit Analyzer(const BLSConfig& config);
  ~Analyzer();

  void setSelection(const std::vector<int>& indices, int natoms);

  bool processFrame(const Frame& frame, FrameMetrics& metrics, std::string& err);

  double gridSpacing() const { return config_.gridSpacing; }

 private:
  BLSConfig config_;
  std::vector<int> selection_;
  bool selectionIsAll_{true};

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace bls
