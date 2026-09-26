#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bls/BLS.hpp"
#include "bls/Comparison.hpp"
#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "config/Parser.hpp"
#include "grid/Grid.hpp"
#include "grid/GridSpec.hpp"
#include "io/Topology.hpp"
#include "io/TrajectoryReader.hpp"
#include "util/Logging.hpp"
#include "util/RSS.hpp"
#include "util/Sha256.hpp"
#include "util/Timer.hpp"

namespace bls {

namespace {

void printUsage() {
  std::cout << "Usage: bls_analyze --system traj.xtc --conf bls.in [options]\n"
               "Options:\n"
               "  --top PATH             Topology file (.gro/.pdb/.xyz)\n"
               "  --out metrics.csv      CSV output path (defaults to stdout)\n"
               "  --json metrics.json    JSON lines output path\n"
               "  --bench bench.csv      Benchmark output path\n"
               "  --stride N             Process every Nth frame\n"
               "  --start N              Skip frames before index N\n"
               "  --stop N               Stop after frame index N (inclusive)\n"
               "  --threads N            Number of OpenMP threads (if enabled)\n"
               "  --replicate N          Replicate ordinal, written to the CSV (default: 1)\n"
               "  --format F             Override molecular system format (xtc,trr,gro,pdb,xyz,mol,sdf)\n"
               "  --algo ALGORITHM       Clustering algorithm to use:\n"
               "                           bls (default) - Bravais Lattice Sampling\n"
               "                           traditional_dfs - Traditional DFS\n"
               "                           skip_dfs - Skip-DFS without lattice\n"
               "                           dbscan - DBSCAN clustering\n"
               "                           hierarchical - Single-linkage hierarchical\n"
               "                           kmeans - K-means clustering\n"
               "                           gcbd - Union-Find grid connectivity\n"
               "                           hdbscan - Hierarchical DBSCAN\n"
               "                           cc3d - Connected Components 3D (fair: basic Union-Find)\n"
               "                           cc3d_optimized - CC3D, SAUF decision-tree two-pass scan\n"
               "                           rle_ccl - Run-Length Encoding CCL (fair: per-voxel union-find)\n"
               "                           rle_ccl_optimized - RLE-CCL with runs as the union-find domain\n"
               "                           vccs - Voxel Cloud Connected Segmentation (fair: uniform seeds)\n"
               "                           vccs_optimized - VCCS with adaptive seeding and seed pruning\n"
               "  --algo-skip N          Jump distance for skip_dfs (default: 3)\n"
               "  --algo-eps F           Epsilon for DBSCAN (default: 3.0)\n"
               "  --algo-minpts N        MinPts for DBSCAN (default: 10)\n"
               "  --algo-k N             K for k-means (default: 20)\n"
               "  --algo-threshold F     Threshold for hierarchical (default: 4.0)\n"
               "  --algo-minclustersize N  Minimum cluster size for HDBSCAN (default: 5)\n"
               "  --algo-minsamples N    Minimum samples for HDBSCAN (default: 5)\n"
               "  --algo-connectivity N  Connectivity for CC3D: 6 or 26 (default: 6)\n"
               "  --compare-plumed PATH  Reference PLUMED CSV/COLVAR\n"
               "  --quiet                Reduce logging\n"
               "  --help                 Show this message\n";
}

double safeParseDouble(const std::string& s, const std::string& opt) {
  try {
    return std::stod(s);
  } catch (...) {
    throw std::runtime_error("Invalid numeric value for " + opt + ": " + s);
  }
}

std::size_t safeParseSize(const std::string& s, const std::string& opt) {
  try {
    return static_cast<std::size_t>(std::stoull(s));
  } catch (...) {
    throw std::runtime_error("Invalid integer value for " + opt + ": " + s);
  }
}

int safeParseInt(const std::string& s, const std::string& opt) {
  try {
    return std::stoi(s);
  } catch (...) {
    throw std::runtime_error("Invalid integer value for " + opt + ": " + s);
  }
}

std::vector<int> buildSelection(const BLSConfig& config, const Topology* topo, int natoms,
                                std::string& err) {
  std::vector<int> indices;
  switch (config.group.type) {
    case GroupSelectorType::All:
      return indices;
    case GroupSelectorType::IndexRange: {
      for (const auto& range : config.group.ranges) {
        int begin = static_cast<int>(range.begin);
        int end = static_cast<int>(range.end);
        for (int idx = begin; idx <= end; ++idx) {
          if (idx >= 0 && idx < natoms) {
            indices.push_back(idx);
          } else {
            Logger::warn("Skipping index ", idx + 1, " outside [1,", natoms, "]");
          }
        }
      }
      return indices;
    }
    case GroupSelectorType::Name: {
      if (!topo) {
        err = "GROUP ATOMS=name requires a topology file.";
        return {};
      }
      for (const auto& name : config.group.names) {
        bool found = false;
        for (const auto& atom : topo->atoms) {
          if (atom.name == name) {
            if (atom.index < natoms) {
              indices.push_back(atom.index);
              found = true;
            }
          }
        }
        if (!found) {
          Logger::warn("No atoms found with name ", name, " in topology.");
        }
      }
      if (indices.empty()) {
        err = "No atoms matched GROUP ATOMS=name selection.";
      }
      return indices;
    }
  }
  return indices;
}

// `replicate` is APPENDED as column 16. The first fifteen columns keep their
// order -- the run scripts read nclusters and elapsed_ms by position (cut -f12,
// -f15), so inserting anywhere but the end would silently shift what those reads
// return. Columns 17-26 were appended by the 2026-09-26 audit: total_ms (the
// whole-frame scope elapsed_ms had before the audit; elapsed_ms is now the labelling
// call alone), the grid fingerprint (h_x..occ_sha256) and BLS's probe count.
void writeCsvHeader(std::ostream& os) {
  os << "frame,time_ps,natoms,NX,NY,NZ,dNN_vox,lattice,centering,seeds,seed_hits,nclusters,"
        "max_cluster,refined_voxels,elapsed_ms,replicate,total_ms,h_x,h_y,h_z,origin_x,"
        "origin_y,origin_z,pbc,occ_sha256,probes\n";
}

void writeCsvRow(std::ostream& os, const FrameMetrics& m, std::size_t frameNumber,
                 int replicate) {
  os << frameNumber << ',' << m.timePs << ',' << m.natoms << ',' << m.nx << ',' << m.ny << ','
     << m.nz << ',' << m.dnnVoxel << ',' << m.lattice << ',' << m.centering << ',' << m.seeds
     << ',' << m.seedHits << ',' << m.nclusters << ',' << m.maxCluster << ','
     << m.refinedVoxels << ',' << m.elapsedMs << ',' << replicate << ',' << m.totalMs << ','
     << m.hx << ',' << m.hy << ',' << m.hz << ',' << m.originX << ',' << m.originY << ','
     << m.originZ << ',' << m.pbc << ',' << m.occSha256 << ',' << m.probes << '\n';
}

void writeJson(std::ostream& os, const FrameMetrics& m, std::size_t frameNumber) {
  os << "{"
     << "\"frame\":" << frameNumber << ","
     << "\"time_ps\":" << m.timePs << ","
     << "\"natoms\":" << m.natoms << ","
     << "\"NX\":" << m.nx << ","
     << "\"NY\":" << m.ny << ","
     << "\"NZ\":" << m.nz << ","
     << "\"dNN_vox\":" << m.dnnVoxel << ","
     << "\"lattice\":\"" << m.lattice << "\","
     << "\"centering\":\"" << m.centering << "\","
     << "\"seeds\":" << m.seeds << ","
     << "\"seed_hits\":" << m.seedHits << ","
     << "\"nclusters\":" << m.nclusters << ","
     << "\"max_cluster\":" << m.maxCluster << ","
     << "\"refined_voxels\":" << m.refinedVoxels << ","
     << "\"elapsed_ms\":" << m.elapsedMs << ","
     << "\"total_ms\":" << m.totalMs << ","
     << "\"h\":[" << m.hx << ',' << m.hy << ',' << m.hz << "],"
     << "\"origin\":[" << m.originX << ',' << m.originY << ',' << m.originZ << "],"
     << "\"pbc\":\"" << m.pbc << "\","
     << "\"occ_sha256\":\"" << m.occSha256 << "\","
     << "\"probes\":" << m.probes;
  if (!m.clusterSizes.empty()) {
    os << ",\"cluster_sizes\":[";
    for (std::size_t i = 0; i < m.clusterSizes.size(); ++i) {
      os << m.clusterSizes[i];
      if (i + 1 < m.clusterSizes.size()) os << ',';
    }
    os << "]";
  }
  os << "}\n";
}

}  // namespace

}  // namespace bls

int main(int argc, char** argv) {
  using namespace bls;

  if (argc == 1) {
    printUsage();
    return EXIT_SUCCESS;
  }

  ProgramOptions opts;
  BLSConfig config;

  auto requireArg = [&](int& i, const std::string& opt) -> std::string {
    if (i + 1 >= argc) {
      throw std::runtime_error("Missing value for option " + opt);
    }
    return std::string(argv[++i]);
  };

  try {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--system") {
        opts.systemPath = requireArg(i, arg);
      } else if (arg == "--top") {
        opts.topologyPath = requireArg(i, arg);
      } else if (arg == "--conf") {
        opts.configPath = requireArg(i, arg);
      } else if (arg == "--out") {
        opts.outputCsvPath = requireArg(i, arg);
      } else if (arg == "--json") {
        opts.outputJsonPath = requireArg(i, arg);
      } else if (arg == "--bench") {
        opts.benchCsvPath = requireArg(i, arg);
      } else if (arg == "--stride") {
        opts.stride = safeParseSize(requireArg(i, arg), arg);
        opts.strideSet = true;
      } else if (arg == "--start") {
        opts.startFrame = safeParseSize(requireArg(i, arg), arg);
      } else if (arg == "--stop") {
        opts.stopFrame = safeParseSize(requireArg(i, arg), arg);
      } else if (arg == "--replicate") {
        opts.replicate = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--threads") {
        opts.threads = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--format") {
        opts.formatOverride = requireArg(i, arg);
      } else if (arg == "--algo") {
        opts.algorithmOverride = requireArg(i, arg);
      } else if (arg == "--algo-skip") {
        opts.skipDfsJumpDistance = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-eps") {
        opts.algoEps = safeParseDouble(requireArg(i, arg), arg);
      } else if (arg == "--algo-minpts") {
        opts.algoMinPts = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-k") {
        opts.algoK = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-threshold") {
        opts.algoThreshold = safeParseDouble(requireArg(i, arg), arg);
      } else if (arg == "--algo-minclustersize") {
        opts.algoMinClusterSize = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-minsamples") {
        opts.algoMinSamples = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-connectivity") {
        opts.algoConnectivity = safeParseInt(requireArg(i, arg), arg);
        opts.algoConnectivitySet = true;
      } else if (arg == "--compare-plumed") {
        opts.comparePlumedPath = requireArg(i, arg);
      } else if (arg == "--quiet") {
        opts.quiet = true;
        Logger::setLevel(LogLevel::Warn);
      } else if (arg == "--help") {
        printUsage();
        return EXIT_SUCCESS;
      } else {
        throw std::runtime_error("Unknown option: " + arg);
      }
    }
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    printUsage();
    return EXIT_FAILURE;
  }

  if (opts.systemPath.empty() || opts.configPath.empty()) {
    std::cerr << "Error: --system and --conf are required.\n";
    return EXIT_FAILURE;
  }

  Parser parser;
  std::string err;
  if (!parser.parseFile(opts.configPath, config, err)) {
    std::cerr << "Config error: " << err << "\n";
    return EXIT_FAILURE;
  }

  std::size_t stride = opts.strideSet ? opts.stride : static_cast<std::size_t>(config.stride);
  if (stride == 0) stride = 1;

  Topology topo;
  Topology* topoPtr = nullptr;
  if (!opts.topologyPath.empty()) {
    if (!loadTopology(opts.topologyPath, topo, err)) {
      std::cerr << "Topology error: " << err << "\n";
      return EXIT_FAILURE;
    }
    topoPtr = &topo;
  }

  auto reader = makeTrajectoryReader(opts.systemPath, opts.formatOverride, err);
  if (!reader) {
    std::cerr << "Molecular system error: " << err << "\n";
    return EXIT_FAILURE;
  }

#ifdef _OPENMP
  if (opts.threads > 0) {
    omp_set_num_threads(opts.threads);
  }
#else
  if (opts.threads > 1) {
    std::cerr << "Warning: binary built without OpenMP support; --threads ignored.\n";
  }
#endif

  std::unique_ptr<std::ofstream> csvFile;
  std::unique_ptr<std::ofstream> jsonFile;
  std::unique_ptr<std::ofstream> benchFile;
  std::ostream* csvStream = &std::cout;
  std::ostream* jsonStream = nullptr;
  std::ostream* benchStream = nullptr;

  if (!opts.outputCsvPath.empty()) {
    csvFile = std::make_unique<std::ofstream>(opts.outputCsvPath);
    if (!*csvFile) {
      std::cerr << "Unable to open CSV output: " << opts.outputCsvPath << "\n";
      return EXIT_FAILURE;
    }
    csvStream = csvFile.get();
  }

  if (!opts.outputJsonPath.empty()) {
    jsonFile = std::make_unique<std::ofstream>(opts.outputJsonPath);
    if (!*jsonFile) {
      std::cerr << "Unable to open JSON output: " << opts.outputJsonPath << "\n";
      return EXIT_FAILURE;
    }
    jsonStream = jsonFile.get();
  }

  if (!opts.benchCsvPath.empty()) {
    benchFile = std::make_unique<std::ofstream>(opts.benchCsvPath);
    if (!*benchFile) {
      std::cerr << "Unable to open bench output: " << opts.benchCsvPath << "\n";
      return EXIT_FAILURE;
    }
    benchStream = benchFile.get();
    *benchStream << "frame,cumulative_ms,peak_rss_bytes,replicate\n";
  }

  writeCsvHeader(*csvStream);

  // Parse and validate algorithm selection
  ClusterAlgorithm selectedAlgo;
  try {
    selectedAlgo = parseAlgorithm(opts.algorithmOverride);
  } catch (const std::exception& ex) {
    std::cerr << "Algorithm error: " << ex.what() << "\n";
    std::cerr << "Available algorithms: ";
    for (const auto& a : listAlgorithms()) std::cerr << a << " ";
    std::cerr << "\n";
    return EXIT_FAILURE;
  }

  if (!opts.quiet) {
    Logger::info("Using clustering algorithm: ", algorithmToString(selectedAlgo));
  }

  Analyzer analyzer(config);
  Grid grid;
  bool selectionReady = false;
  std::vector<int> selection;
  std::vector<FrameMetrics> frames;

  std::size_t frameIndex = 0;
  double cumulativeMs = 0.0;
  std::size_t peakRss = 0;

  while (true) {
    Frame current;
    if (!reader->read(current, err)) {
      if (!err.empty()) {
        std::cerr << "System read error: " << err << "\n";
        return EXIT_FAILURE;
      }
      break;
    }

    int natoms = current.natoms;
    if (!selectionReady) {
      std::string selErr;
      selection = buildSelection(config, topoPtr, natoms, selErr);
      if (!selErr.empty()) {
        std::cerr << "Selection error: " << selErr << "\n";
        return EXIT_FAILURE;
      }
      try {
        analyzer.setSelection(selection, natoms);
      } catch (const std::exception& ex) {
        std::cerr << "Selection error: " << ex.what() << "\n";
        return EXIT_FAILURE;
      }
      selectionReady = true;
    }

    if (frameIndex < opts.startFrame) {
      ++frameIndex;
      continue;
    }
    if (frameIndex > opts.stopFrame) {
      break;
    }
    if ((frameIndex - opts.startFrame) % stride != 0) {
      ++frameIndex;
      continue;
    }

    FrameMetrics metrics;

    // One grid for every method (Defect 12 fix): the box, origin, dimensions and
    // periodicity come from the shared builder; the voxelisation is identical; the
    // labelling timer starts only after both, for BLS and for every comparison method.
    ScopedTimer totalTimer;
    GridSpec spec;
    const std::vector<int>* selPtr = selection.empty() ? nullptr : &selection;
    if (!deriveGridSpec(config, current, selPtr, spec, err)) {
      std::cerr << "Error: " << err << "\n";
      return EXIT_FAILURE;
    }
    configureGrid(grid, spec);
    grid.rasterize(current.xyz, selPtr, config.cutoff, config.occupancy);

    if (selectedAlgo == ClusterAlgorithm::BLS) {
      if (!analyzer.labelGrid(grid, metrics, err)) {
        std::cerr << "Processing error: " << err << "\n";
        return EXIT_FAILURE;
      }
    } else {
      if (spec.pbc.any()) {
        std::cerr << "Error: PBC is not implemented for " << algorithmToString(selectedAlgo)
                  << "\n";
        return EXIT_FAILURE;
      }
      ClusterParams params;
      params.nx = spec.nx;
      params.ny = spec.ny;
      params.nz = spec.nz;
      params.skipDfsJumpDistance = opts.skipDfsJumpDistance;
      params.eps = opts.algoEps;
      params.minPts = opts.algoMinPts;
      params.k = opts.algoK;
      params.threshold = opts.algoThreshold;
      params.connectivity =
          opts.algoConnectivitySet ? opts.algoConnectivity : config.connectivity;
      params.minClusterSize = opts.algoMinClusterSize;
      params.minSamples = opts.algoMinSamples;

      ScopedTimer labelTimer;
      ClusterResult result =
          runClusterAlgorithm(selectedAlgo, params, grid.occupancy(), grid.visited());
      metrics.elapsedMs = labelTimer.elapsedMilliseconds();

      metrics.dnnVoxel = 0.0;  // Not applicable for non-BLS algorithms
      metrics.lattice = algorithmToString(selectedAlgo);
      metrics.centering = "-";
      metrics.seeds = 0;
      metrics.seedHits = 0;
      metrics.nclusters = result.nclusters;
      metrics.maxCluster = result.maxCluster;
      metrics.refinedVoxels = result.visitedVoxels;
      metrics.clusterSizes = std::move(result.clusterSizes);
    }
    metrics.totalMs = totalTimer.elapsedMilliseconds();

    metrics.timePs = current.time;
    metrics.natoms = current.natoms;
    metrics.nx = spec.nx;
    metrics.ny = spec.ny;
    metrics.nz = spec.nz;
    // Fingerprint, outside both timers.
    metrics.hx = spec.h(0);
    metrics.hy = spec.h(1);
    metrics.hz = spec.h(2);
    metrics.originX = spec.origin.x;
    metrics.originY = spec.origin.y;
    metrics.originZ = spec.origin.z;
    metrics.pbc = spec.pbc.toString();
    metrics.occSha256 = occupancySha256(grid.occupancy());

    metrics.frameIndex = frameIndex;

    writeCsvRow(*csvStream, metrics, frameIndex, opts.replicate);
    if (jsonStream) {
      writeJson(*jsonStream, metrics, frameIndex);
    }

    cumulativeMs += metrics.elapsedMs;
    peakRss = std::max(peakRss, currentRSSBytes());
    if (benchStream) {
      *benchStream << frameIndex << ',' << cumulativeMs << ',' << peakRss << ','
                   << opts.replicate << '\n';
    }

    frames.push_back(metrics);
    ++frameIndex;
  }

  reader->close();

  if (!opts.comparePlumedPath.empty()) {
    ComparisonSummary summary;
    if (compareWithPlumed(opts.comparePlumedPath, frames, summary, err)) {
      std::cout << "# PLUMED comparison over " << summary.matchedFrames << " frames\n"
                << "# mean|max_cluster| diff: " << summary.meanAbsMaxCluster
                << ", rmse: " << summary.rmseMaxCluster << '\n'
                << "# mean|nclusters| diff: " << summary.meanAbsNClusters
                << ", rmse: " << summary.rmseNClusters << '\n'
                << "# mean elapsed diff (ms): " << summary.meanElapsedDiff
                << ", speedup: " << summary.speedup << "x\n"
                << "# Kendall tau (cluster sizes): " << summary.kendallTau << '\n';
    } else {
      std::cerr << "Comparison error: " << err << "\n";
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
