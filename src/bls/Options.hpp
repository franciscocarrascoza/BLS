#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace bls {

struct ProgramOptions {
  std::string systemPath;
  std::string topologyPath;
  std::string configPath;
  std::string outputCsvPath;
  std::string outputJsonPath;
  std::string benchCsvPath;
  std::string comparePlumedPath;
  std::string formatOverride{"auto"};
  std::string algorithmOverride{"bls"};  // Clustering algorithm selection
  std::size_t stride{1};
  bool strideSet{false};
  std::size_t startFrame{0};
  std::size_t stopFrame{std::numeric_limits<std::size_t>::max()};
  int threads{1};
  bool wantBench{false};
  bool quiet{false};
  // Algorithm-specific parameters
  int algoSkip{3};           // Skip distance for skip_dfs
  double algoEps{3.0};       // Epsilon for DBSCAN
  int algoMinPts{10};        // MinPts for DBSCAN
  int algoK{20};             // K for k-means/spectral
  double algoThreshold{4.0}; // Threshold for hierarchical
  int algoMinClusterSize{5}; // Minimum cluster size for HDBSCAN
  int algoMinSamples{5};     // Minimum samples for HDBSCAN
  int algoConnectivity{6};   // Connectivity for CC3D (6 or 26)
};

enum class GroupSelectorType { All, IndexRange, Name };

struct IndexRange {
  std::size_t begin{0};
  std::size_t end{0};  // inclusive
};

struct AtomSelection {
  GroupSelectorType type{GroupSelectorType::All};
  std::vector<IndexRange> ranges;
  std::vector<std::string> names;
};

enum class BoxMode { Auto, Manual };

// Whether the analysis box is to be treated as a periodic cell.
//
// This has to be explicit because the two halves of Grid::rasterize used to
// disagree about it: Grid::fractional wrapped atom positions with floor()
// (periodic) while the neighbour stencil clipped out-of-range voxel indices
// (non-periodic). An atom near a face therefore had its position folded into
// the cell but only half its footprint stamped. Both halves now read this
// flag, so a box is periodic in both or neither.
//
// NonPeriodic is right for a box the analyser synthesised itself: BLS.cpp
// builds a bounding box around the atoms with 2*GRID_SPACING of padding when
// the trajectory carries no usable cell, and wrapping that would connect
// opposite faces across vacuum. Periodic is right for a genuine MD cell.
enum class BoxPeriodicity { NonPeriodic, Periodic };

// Single place where the periodicity of the analysis box is decided, used by
// both box-derivation sites (bls/BLS.cpp and main.cpp).
//
// CAVEAT, deliberate and worth knowing before changing this: BOX AUTO does not
// always mean a synthesised bounding box. It uses the trajectory's own cell
// whenever that cell is usable, and only falls back to a bounding box when the
// cell is missing, degenerate, or implausibly large. So a genuine periodic MD
// cell arriving through a CRYST1 record is treated as non-periodic here. That
// costs nothing today -- no clustering algorithm in this codebase has any
// periodic-boundary handling at all, and in every E0-E5 system every atom lies
// strictly inside its cell, so no wrap would fire even if it were enabled --
// but it is the line to revisit if periodic connectivity is ever added.
inline BoxPeriodicity periodicityForBoxMode(BoxMode mode) {
  return mode == BoxMode::Manual ? BoxPeriodicity::Periodic : BoxPeriodicity::NonPeriodic;
}

struct ManualBox {
  double xlo{0.0}, xhi{0.0};
  double ylo{0.0}, yhi{0.0};
  double zlo{0.0}, zhi{0.0};
};

enum class OccupancyMode { Any, All };

enum class LatticeType { Cubic, Hexagonal, Triclinic };
enum class CenteringType { P, F, I };

struct LatticeSettings {
  LatticeType lattice{LatticeType::Cubic};
  CenteringType centering{CenteringType::F};
  double hexCOverA{1.633};
  double triclinicA{1.0};
  double triclinicB{1.2};
  double triclinicC{1.4};
  double triclinicAlphaDeg{90.0};
  double triclinicBetaDeg{100.0};
  double triclinicGammaDeg{110.0};
};

struct BLSConfig {
  AtomSelection group;
  BoxMode boxMode{BoxMode::Auto};
  ManualBox manualBox{};
  double gridSpacing{0.25};
  int connectivity{6};
  int skip{3};
  double alpha{0.7};
  double dnn{0.0};
  bool hasExplicitDnn{false};
  std::vector<double> radii;
  double cutoff{0.0};
  OccupancyMode occupancy{OccupancyMode::Any};
  LatticeSettings lattice;
  int stride{1};
  std::vector<std::string> outputs{
      "NCLUSTERS", "MAX_CLUSTER", "SEED_HITS", "SEEDS", "REFINED_VOXELS"};
};

struct InputDeck {
  ProgramOptions cli;
  BLSConfig config;
};

}  // namespace bls
