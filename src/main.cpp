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
#include "config/Parser.hpp"
#include "io/Topology.hpp"
#include "io/TrajectoryReader.hpp"
#include "util/Logging.hpp"
#include "util/RSS.hpp"

namespace bls {

namespace {

void printUsage() {
  std::cout << "Usage: bls_analyze --traj traj.xtc --conf bls.in [options]\n"
               "Options:\n"
               "  --top PATH             Topology file (.gro/.pdb)\n"
               "  --out metrics.csv      CSV output path (defaults to stdout)\n"
               "  --json metrics.json    JSON lines output path\n"
               "  --bench bench.csv      Benchmark output path\n"
               "  --stride N             Process every Nth frame\n"
               "  --start N              Skip frames before index N\n"
               "  --stop N               Stop after frame index N (inclusive)\n"
               "  --threads N            Number of OpenMP threads (if enabled)\n"
               "  --format F             Override trajectory format\n"
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
    if (!analyzer.processFrame(current, metrics, err)) {
      std::cerr << "Processing error: " << err << "\n";
      return EXIT_FAILURE;
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
