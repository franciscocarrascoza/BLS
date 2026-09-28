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
  // Replicate ordinal, emitted as the last CSV column so the campaign's repeats
  // land in one file. Purely a label: it changes no computation.
  int replicate{1};
  // Algorithm-specific parameters.
  //
  // One name per meaning, deliberately. Until Task 9 a single --algo-skip fed
  // two unrelated algorithms: skip_dfs's jump distance and octree_ccl's leaf
  // size, which is how octree_ccl silently ran at leaf 3 against its own
  // documented intent of 8. octree_ccl was removed in Task 12; the naming rule
  // it motivated stays. See also BLSConfig::refinementStride, a third distinct
  // quantity that was also called "skip".
  int skipDfsJumpDistance{3};  // --algo-skip: jump distance for skip_dfs
  double algoEps{3.0};       // Epsilon for DBSCAN
  int algoMinPts{10};        // MinPts for DBSCAN
  int algoK{20};             // K for k-means
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

// AUTO   -- the trajectory's own cell when usable, else a bounding box fitted to the
//           (selected) atoms with 2*GRID_SPACING padding; see deriveGridSpec().
// MANUAL -- the orthorhombic box given by XLO..ZHI in the deck.
// CELL   -- the frame's own cell (PDB CRYST1, per MODEL), origin (0,0,0); an error if the
//           frame carries no usable cell. Required for periodic runs on PDB input.
// Every method -- BLS and all comparison algorithms -- grids the box chosen here through
// the one shared builder, deriveGridSpec() (grid/GridSpec.hpp). Before the 2026-09-26
// audit, main.cpp carried its own copy of the AUTO logic and never read MANUAL
// (Defect 12), so BLS and the comparison methods gridded different boxes.
enum class BoxMode { Auto, Manual, Cell };

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

// Periodicity is stated by the deck (keyword PBC), never inferred from the box mode.
//   PBC none  (default) -- non-periodic in every method; reproduces pre-audit outputs.
//   PBC xyz             -- periodic along all three cell vectors, in the voxeliser and in
//                          every method's adjacency. Requires BOX CELL or BOX MANUAL: a
//                          bounding box synthesised around the atoms is not a periodic cell.
// Per-axis periodicity (e.g. PBC xy) is rejected by the parser: no method implements it.
//
// History: until the audit, BOX MANUAL silently implied periodic RASTERISATION while
// no clustering algorithm had any periodic adjacency, and BOX AUTO treated a genuine
// CRYST1 cell as non-periodic. Both halves now read this one flag.
struct PeriodicAxes {
  bool x{false}, y{false}, z{false};
  bool any() const { return x || y || z; }
  bool all() const { return x && y && z; }
  std::string toString() const {
    if (!any()) return "none";
    std::string s;
    if (x) s += 'x';
    if (y) s += 'y';
    if (z) s += 'z';
    return s;
  }
};

inline BoxPeriodicity periodicityFor(const PeriodicAxes& pbc) {
  return pbc.all() ? BoxPeriodicity::Periodic : BoxPeriodicity::NonPeriodic;
}

struct ManualBox {
  double xlo{0.0}, xhi{0.0};
  double ylo{0.0}, yhi{0.0};
  double zlo{0.0}, zhi{0.0};
};

enum class OccupancyMode { Any, All };

// BLS's refinement (stage 2, after the probe evaluation). Deck keyword REFINEMENT.
//   SKIP_DFS (default) -- refine::SkipDFS, advancing SKIP voxels per step along each direction;
//                         the method as the manuscript describes it. At SKIP 1 it is a plain
//                         DFS that marks voxels when they are pushed.
//   DFS                -- refine::StandardDFS, the traditional_dfs baseline's own flood fill
//                         (mark when popped), started from the seeds instead of from a raster
//                         scan. Added by the audit (28.09.26) so that BLS - DFS isolates the
//                         seeding: with DFS refinement the two methods differ only in how the
//                         components are found (probe sites vs every voxel).
// Both return the same components, sizes and labels (6-connectivity: a skip step stops at the
// first unoccupied voxel, so it only reaches voxels joined through occupied ones).
enum class RefinementMode { SkipDFS, DFS };

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
  PeriodicAxes pbc{};  // deck keyword PBC; default none
  // Deck keyword LATTICE_ORIGIN ox oy oz: translation of the BLS probe lattice, in
  // voxels, applied to every site before it is rounded to a voxel (modulo the cell
  // when periodic). Used to sweep lattice-origin offsets without moving atoms, which
  // would change the rasterised occupancy. Default (0,0,0) = the historical anchoring.
  double latticeOrigin[3]{0.0, 0.0, 0.0};
  double gridSpacing{0.25};
  int connectivity{6};
  // Config keyword SKIP. BLS's refinement (refine::SkipDFS) advances this many
  // voxels per step; it is NOT ProgramOptions::skipDfsJumpDistance, which
  // belongs to the unrelated cluster::skipDFS comparison algorithm.
  int refinementStride{3};
  RefinementMode refinement{RefinementMode::SkipDFS};  // deck keyword REFINEMENT
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
