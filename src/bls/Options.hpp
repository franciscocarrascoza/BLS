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
  bool quiet{false};
  // Algorithm-specific parameters.
  //
  // One name per meaning, deliberately. Until Task 9 a single --algo-skip fed
  // two unrelated algorithms: skip_dfs's jump distance and octree_ccl's leaf
  // size. octree_ccl's own comment documented an intended leaf size of 8, but
  // it inherited algoSkip's default of 3 instead, so every published
  // octree_ccl timing ran a ~6% slower tree than the one described. They are
  // separate parameters now; see also BLSConfig::refinementStride, a third
  // distinct quantity that was also called "skip".
  int skipDfsJumpDistance{3};  // --algo-skip: jump distance for skip_dfs
  int octreeLeafSize{8};       // --octree-leaf: leaf edge, voxels, octree_ccl
  double algoEps{3.0};       // Epsilon for DBSCAN
  int algoMinPts{10};        // MinPts for DBSCAN
  int algoK{20};             // K for k-means/spectral
  double algoThreshold{4.0}; // Threshold for hierarchical
  int algoMinClusterSize{5}; // Minimum cluster size for HDBSCAN
  int algoMinSamples{5};     // Minimum samples for HDBSCAN
  // 6 is a legal value, so it cannot also mean "unset" -- with the old
  // sentinel test (algoConnectivity != 6) an explicit --algo-connectivity 6
  // was indistinguishable from no flag at all, and silently lost to a config
  // saying 26. Tracked with an explicit flag, as --stride already does.
  int algoConnectivity{6};   // Connectivity for CC3D (6 or 26)
  bool algoConnectivitySet{false};
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
  // No default worth trusting. cubic/F used to be the silent fallback, and
  // while the enumerator was broken (fixed in a90066e) that fallback was inert
  // -- every occupied voxel became a seed regardless of lattice, so the field
  // did nothing. It is load-bearing now: on ld-asw, cubic/F versus the
  // structure's own lattice moves the cluster count by 17%. A config that
  // omits LATTICE or CENTERING is therefore rejected rather than defaulted;
  // `set` records whether the config said so. The initialisers below are only
  // the pre-parse state and must never be read as a policy choice.
  LatticeType lattice{LatticeType::Cubic};
  bool latticeSet{false};
  CenteringType centering{CenteringType::F};
  bool centeringSet{false};
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
  // Config keyword SKIP. BLS's refinement (refine::SkipDFS) advances this many
  // voxels per step; it is NOT ProgramOptions::skipDfsJumpDistance, which
  // belongs to the unrelated cluster::skipDFS comparison algorithm.
  int refinementStride{3};
  double alpha{0.7};
  double dnn{0.0};
  bool hasExplicitDnn{false};
  std::vector<double> radii;
  double cutoff{0.0};
  OccupancyMode occupancy{OccupancyMode::Any};
  LatticeSettings lattice;
  int stride{1};
  // No `outputs` field: the CSV and JSON column sets are fixed in main.cpp.
  // The OUTPUT keyword that used to fill this was parsed and then never read
  // by anything; see Parser.cpp for why it was deleted rather than wired up.
};

struct InputDeck {
  ProgramOptions cli;
  BLSConfig config;
};

}  // namespace bls
