#pragma once

#include <string>
#include <vector>

#include "bls/BLS.hpp"

namespace bls {

struct ComparisonSummary {
  double meanAbsMaxCluster{0.0};
  double rmseMaxCluster{0.0};
  double meanAbsNClusters{0.0};
  double rmseNClusters{0.0};
  double meanElapsedDiff{0.0};
  double speedup{0.0};
  double kendallTau{1.0};
  std::size_t matchedFrames{0};
};

bool compareWithPlumed(const std::string& path, const std::vector<FrameMetrics>& ours,
                       ComparisonSummary& summary, std::string& err);

}  // namespace bls

