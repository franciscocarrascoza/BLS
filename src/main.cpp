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
#include "io/Topology.hpp"
#include "io/TrajectoryReader.hpp"
#include "util/Logging.hpp"
#include "util/RSS.hpp"
#include "util/Timer.hpp"

namespace bls {

namespace {

void printUsage() {
  std::cout << "Usage: bls_analyze --traj traj.xtc --conf bls.in [options]\n"
               "Options:\n"
               "  --top PATH             Topology file (.gro/.pdb/.xyz)\n"
               "  --out metrics.csv      CSV output path (defaults to stdout)\n"
               "  --json metrics.json    JSON lines output path\n"
               "  --bench bench.csv      Benchmark output path\n"
               "  --stride N             Process every Nth frame\n"
               "  --start N              Skip frames before index N\n"
               "  --stop N               Stop after frame index N (inclusive)\n"
               "  --threads N            Number of OpenMP threads (if enabled)\n"
               "  --format F             Override trajectory format (xtc,trr,gro,pdb,xyz,mol,sdf)\n"
               "  --algo ALGORITHM       Clustering algorithm to use:\n"
               "                           bls (default) - Bravais Lattice Sampling\n"
               "                           traditional_dfs - Traditional DFS\n"
               "                           skip_dfs - Skip-DFS without lattice\n"
               "                           dbscan - DBSCAN clustering\n"
               "                           hierarchical - Single-linkage hierarchical\n"
               "                           kmeans - K-means clustering\n"
               "                           spectral - Spectral clustering (simplified)\n"
               "                           gcbd - Union-Find grid connectivity\n"
               "  --algo-skip N          Skip distance for skip_dfs (default: 3)\n"
               "  --algo-eps F           Epsilon for DBSCAN (default: 3.0)\n"
               "  --algo-minpts N        MinPts for DBSCAN (default: 10)\n"
               "  --algo-k N             K for k-means/spectral (default: 20)\n"
               "  --algo-threshold F     Threshold for hierarchical (default: 4.0)\n"
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

void writeCsvHeader(std::ostream& os) {
  os << "frame,time_ps,natoms,NX,NY,NZ,dNN_vox,lattice,centering,seeds,seed_hits,nclusters,"
        "max_cluster,refined_voxels,elapsed_ms\n";
}

void writeCsvRow(std::ostream& os, const FrameMetrics& m, std::size_t frameNumber) {
  os << frameNumber << ',' << m.timePs << ',' << m.natoms << ',' << m.nx << ',' << m.ny << ','
     << m.nz << ',' << m.dnnVoxel << ',' << m.lattice << ',' << m.centering << ',' << m.seeds
     << ',' << m.seedHits << ',' << m.nclusters << ',' << m.maxCluster << ','
     << m.refinedVoxels << ',' << m.elapsedMs << '\n';
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
     << "\"elapsed_ms\":" << m.elapsedMs;
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
      if (arg == "--traj") {
        opts.trajectoryPath = requireArg(i, arg);
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
      } else if (arg == "--threads") {
        opts.threads = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--format") {
        opts.formatOverride = requireArg(i, arg);
      } else if (arg == "--algo") {
        opts.algorithmOverride = requireArg(i, arg);
      } else if (arg == "--algo-skip") {
        opts.algoSkip = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-eps") {
        opts.algoEps = safeParseDouble(requireArg(i, arg), arg);
      } else if (arg == "--algo-minpts") {
        opts.algoMinPts = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-k") {
        opts.algoK = safeParseInt(requireArg(i, arg), arg);
      } else if (arg == "--algo-threshold") {
        opts.algoThreshold = safeParseDouble(requireArg(i, arg), arg);
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

  if (opts.trajectoryPath.empty() || opts.configPath.empty()) {
    std::cerr << "Error: --traj and --conf are required.\n";
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

  auto reader = makeTrajectoryReader(opts.trajectoryPath, opts.formatOverride, err);
  if (!reader) {
    std::cerr << "Trajectory error: " << err << "\n";
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
    *benchStream << "frame,cumulative_ms,peak_rss_bytes\n";
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
        std::cerr << "Trajectory read error: " << err << "\n";
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

    if (selectedAlgo == ClusterAlgorithm::BLS) {
      // Use standard BLS processing via Analyzer
      if (!analyzer.processFrame(current, metrics, err)) {
        std::cerr << "Processing error: " << err << "\n";
        return EXIT_FAILURE;
      }
    } else {
      // Use alternative clustering algorithm
      ScopedTimer timer;

      // Set up grid similar to Analyzer
      Mat3 activeBox = current.box;
      Vec3 origin{0.0, 0.0, 0.0};

      auto col0 = activeBox.column(0);
      auto col1 = activeBox.column(1);
      auto col2 = activeBox.column(2);
      double len0 = norm(col0);
      double len1 = norm(col1);
      double len2 = norm(col2);

      // Helper lambda to compute coordinate bounds
      auto computeBounds = [&](Vec3& minPos, Vec3& maxPos) {
        minPos = Vec3{std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity()};
        maxPos = Vec3{-std::numeric_limits<double>::infinity(),
                      -std::numeric_limits<double>::infinity(),
                      -std::numeric_limits<double>::infinity()};
        for (const auto& p : current.xyz) {
          minPos.x = std::min(minPos.x, p.x);
          minPos.y = std::min(minPos.y, p.y);
          minPos.z = std::min(minPos.z, p.z);
          maxPos.x = std::max(maxPos.x, p.x);
          maxPos.y = std::max(maxPos.y, p.y);
          maxPos.z = std::max(maxPos.z, p.z);
        }
      };

      // Check if box needs correction (zero/invalid or unreasonably large)
      bool needsBoxCorrection = (len0 < 1e-8 || len1 < 1e-8 || len2 < 1e-8);

      if (!needsBoxCorrection && !current.xyz.empty()) {
        Vec3 minPos, maxPos;
        computeBounds(minPos, maxPos);

        double coordExtentX = maxPos.x - minPos.x;
        double coordExtentY = maxPos.y - minPos.y;
        double coordExtentZ = maxPos.z - minPos.z;

        // If box is more than 10x larger than coordinate extent, it's likely incorrect
        const double suspiciousRatio = 10.0;
        bool boxTooLarge = (len0 > coordExtentX * suspiciousRatio) ||
                           (len1 > coordExtentY * suspiciousRatio) ||
                           (len2 > coordExtentZ * suspiciousRatio);

        if (boxTooLarge) {
          Logger::warn("Box size (", len0, " x ", len1, " x ", len2,
                       ") is unreasonably large compared to coordinate extent (",
                       coordExtentX, " x ", coordExtentY, " x ", coordExtentZ,
                       "). Auto-correcting to fit coordinates.");
          needsBoxCorrection = true;
        }
      }

      if (needsBoxCorrection && !current.xyz.empty()) {
        Vec3 minPos, maxPos;
        computeBounds(minPos, maxPos);
        double padding = config.gridSpacing * 2.0;
        origin = minPos;
        activeBox = Mat3{Vec3{std::max(maxPos.x - minPos.x + padding, padding), 0.0, 0.0},
                         Vec3{0.0, std::max(maxPos.y - minPos.y + padding, padding), 0.0},
                         Vec3{0.0, 0.0, std::max(maxPos.z - minPos.z + padding, padding)}};
        col0 = activeBox.column(0);
        col1 = activeBox.column(1);
        col2 = activeBox.column(2);
      }

      int nx = std::max(1, static_cast<int>(std::ceil(norm(col0) / config.gridSpacing)));
      int ny = std::max(1, static_cast<int>(std::ceil(norm(col1) / config.gridSpacing)));
      int nz = std::max(1, static_cast<int>(std::ceil(norm(col2) / config.gridSpacing)));

      // Check memory requirements before allocation
      std::size_t requiredMemory = estimateGridMemoryBytes(nx, ny, nz);
      std::size_t availableMemory = availableSystemRAMBytes();
      std::size_t maxAllowedMemory = static_cast<std::size_t>(availableMemory * 0.8);

      if (requiredMemory > maxAllowedMemory) {
        std::cerr << "Error: Grid allocation would require "
                  << (requiredMemory / (1024.0 * 1024.0 * 1024.0)) << " GB, "
                  << "which exceeds available RAM limit ("
                  << (maxAllowedMemory / (1024.0 * 1024.0 * 1024.0)) << " GB).\n";
        std::cerr << "Grid dimensions: " << nx << " x " << ny << " x " << nz << " = "
                  << (static_cast<std::size_t>(nx) * ny * nz) << " voxels\n";
        std::cerr << "Box size: " << norm(col0) << " x " << norm(col1) << " x " << norm(col2)
                  << " Angstroms\n";
        std::cerr << "Grid spacing: " << config.gridSpacing << " Angstroms\n\n";
        std::cerr << "Solutions:\n";
        std::cerr << "  1. Increase GRID_SPACING (current: " << config.gridSpacing << " A)\n";
        double minSpacing = std::max({norm(col0), norm(col1), norm(col2)}) /
                            maxBoxDimensionForRAM(1.0, maxAllowedMemory);
        std::cerr << "     Minimum spacing for this box: " << minSpacing << " A\n";
        std::cerr << "  2. Reduce box size (check CRYST1 record in PDB or use BOX MANUAL)\n";
        double maxBoxSize = maxBoxDimensionForRAM(config.gridSpacing, maxAllowedMemory);
        std::cerr << "     Maximum box dimension for current spacing: " << maxBoxSize << " A\n";
        return EXIT_FAILURE;
      }

      Grid grid;
      grid.configure(nx, ny, nz, config.gridSpacing, activeBox, origin);
      const std::vector<int>* selPtr = selection.empty() ? nullptr : &selection;
      grid.rasterize(current.xyz, selPtr, config.cutoff, config.occupancy);

      // Set up algorithm parameters
      ClusterParams params;
      params.nx = nx;
      params.ny = ny;
      params.nz = nz;
      params.skip = opts.algoSkip;
      params.eps = opts.algoEps;
      params.minPts = opts.algoMinPts;
      params.k = opts.algoK;
      params.threshold = opts.algoThreshold;
      params.connectivity = config.connectivity;

      // Run the selected algorithm
      ClusterResult result = runClusterAlgorithm(
          selectedAlgo, params, grid.occupancy(), grid.visited());

      // Fill metrics
      metrics.timePs = current.time;
      metrics.natoms = current.natoms;
      metrics.nx = nx;
      metrics.ny = ny;
      metrics.nz = nz;
      metrics.dnnVoxel = 0.0;  // Not applicable for non-BLS algorithms
      metrics.lattice = algorithmToString(selectedAlgo);
      metrics.centering = "-";
      metrics.seeds = 0;
      metrics.seedHits = 0;
      metrics.nclusters = result.nclusters;
      metrics.maxCluster = result.maxCluster;
      metrics.refinedVoxels = result.visitedVoxels;
      metrics.clusterSizes = std::move(result.clusterSizes);
      metrics.elapsedMs = timer.elapsedMilliseconds();
    }

    metrics.frameIndex = frameIndex;

    writeCsvRow(*csvStream, metrics, frameIndex);
    if (jsonStream) {
      writeJson(*jsonStream, metrics, frameIndex);
    }

    cumulativeMs += metrics.elapsedMs;
    peakRss = std::max(peakRss, currentRSSBytes());
    if (benchStream) {
      *benchStream << frameIndex << ',' << cumulativeMs << ',' << peakRss << '\n';
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
