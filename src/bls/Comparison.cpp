#include "bls/Comparison.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace bls {

namespace {

std::string trim(const std::string& s) {
  auto begin = std::find_if_not(s.begin(), s.end(),
                                [](unsigned char c) { return std::isspace(c); });
  auto end = std::find_if_not(s.rbegin(), s.rend(),
                              [](unsigned char c) { return std::isspace(c); })
                 .base();
  if (begin >= end) return {};
  return std::string(begin, end);
}

std::vector<std::string> split(const std::string& line, char delim) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream iss(line);
  if (delim == ' ') {
    while (iss >> token) {
      tokens.push_back(token);
    }
  } else {
    while (std::getline(iss, token, delim)) {
      tokens.push_back(trim(token));
    }
  }
  return tokens;
}

struct ReferenceFrame {
  std::size_t frame{0};
  double time{0.0};
  int nclusters{0};
  int maxCluster{0};
  double elapsedMs{0.0};
  std::vector<int> clusterSizes;
};

struct Columns {
  int frame{-1};
  int time{-1};
  int nclusters{-1};
  int maxCluster{-1};
  int elapsed{-1};
  int clusterSizes{-1};
};

int findColumn(const std::vector<std::string>& headers, const std::string& name) {
  for (std::size_t i = 0; i < headers.size(); ++i) {
    std::string lower = headers[i];
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == name) return static_cast<int>(i);
  }
  return -1;
}

std::vector<int> parseClusterSizes(const std::string& field) {
  std::vector<int> values;
  auto tokens = split(field, ';');
  if (tokens.size() == 1) {
    tokens = split(field, ',');
  }
  for (const auto& tok : tokens) {
    if (tok.empty()) continue;
    values.push_back(std::stoi(tok));
  }
  std::sort(values.begin(), values.end(), std::greater<int>());
  return values;
}

bool loadReference(const std::string& path, std::vector<ReferenceFrame>& out, std::string& err) {
  std::ifstream in(path);
  if (!in) {
    err = "Unable to open PLUMED comparison file: " + path;
    return false;
  }

  Columns cols;
  bool headerParsed = false;
  std::size_t inferredFrame = 0;
  char delim = ',';

  std::string line;
  while (std::getline(in, line)) {
    auto trimmed = trim(line);
    if (trimmed.empty() || trimmed[0] == '@') continue;

    if (!headerParsed) {
      std::string header = trimmed;
      if (trimmed.rfind("#!", 0) == 0) {
        header = trim(trimmed.substr(2));
        if (header.rfind("FIELDS", 0) == 0) {
          header = trim(header.substr(6));
        }
      } else if (trimmed[0] == '#') {
        header = trim(trimmed.substr(1));
      }
      if (!header.empty()) {
        delim = header.find(',') != std::string::npos ? ',' : ' ';
        auto headers = split(header, delim);
        cols.frame = findColumn(headers, "frame");
        cols.time = findColumn(headers, "time");
        if (cols.time < 0) cols.time = findColumn(headers, "t");
        cols.nclusters = findColumn(headers, "nclusters");
        if (cols.nclusters < 0) cols.nclusters = findColumn(headers, "ncluster");
        cols.maxCluster = findColumn(headers, "max_cluster");
        cols.elapsed = findColumn(headers, "elapsed_ms");
        if (cols.elapsed < 0) cols.elapsed = findColumn(headers, "elapsed");
        cols.clusterSizes = findColumn(headers, "cluster_sizes");
        headerParsed = (cols.time >= 0 || cols.nclusters >= 0 || cols.maxCluster >= 0);
        if (headerParsed) continue;
      }
    }

    delim = trimmed.find(',') != std::string::npos ? ',' : delim;
    auto tokens = split(trimmed, delim);
    if (tokens.empty()) continue;

    ReferenceFrame ref;
    if (cols.frame >= 0 && cols.frame < static_cast<int>(tokens.size())) {
      ref.frame = static_cast<std::size_t>(std::stoul(tokens[cols.frame]));
    } else {
      ref.frame = inferredFrame;
    }
    if (cols.time >= 0 && cols.time < static_cast<int>(tokens.size())) {
      ref.time = std::stod(tokens[cols.time]);
    }
    if (cols.nclusters >= 0 && cols.nclusters < static_cast<int>(tokens.size())) {
      ref.nclusters = std::stoi(tokens[cols.nclusters]);
    }
    if (cols.maxCluster >= 0 && cols.maxCluster < static_cast<int>(tokens.size())) {
      ref.maxCluster = std::stoi(tokens[cols.maxCluster]);
    }
    if (cols.elapsed >= 0 && cols.elapsed < static_cast<int>(tokens.size())) {
      ref.elapsedMs = std::stod(tokens[cols.elapsed]);
    }
    if (cols.clusterSizes >= 0 && cols.clusterSizes < static_cast<int>(tokens.size())) {
      ref.clusterSizes = parseClusterSizes(tokens[cols.clusterSizes]);
    }

    out.push_back(ref);
    ++inferredFrame;
  }

  if (out.empty()) {
    err = "No data rows parsed from PLUMED file.";
    return false;
  }
  return true;
}

double kendallTau(const std::vector<int>& a, const std::vector<int>& b) {
  std::size_t n = std::min(a.size(), b.size());
  if (n < 2) return 1.0;
  int concordant = 0;
  int discordant = 0;
  for (std::size_t i = 0; i < n - 1; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      int sa = (a[i] < a[j]) ? -1 : (a[i] > a[j]) ? 1 : 0;
      int sb = (b[i] < b[j]) ? -1 : (b[i] > b[j]) ? 1 : 0;
      if (sa == 0 || sb == 0) continue;
      if (sa == sb)
        ++concordant;
      else
        ++discordant;
    }
  }
  int denom = concordant + discordant;
  if (denom == 0) return 1.0;
  return static_cast<double>(concordant - discordant) / static_cast<double>(denom);
}

}  // namespace

bool compareWithPlumed(const std::string& path, const std::vector<FrameMetrics>& ours,
                       ComparisonSummary& summary, std::string& err) {
  std::vector<ReferenceFrame> ref;
  if (!loadReference(path, ref, err)) {
    return false;
  }

  std::size_t count = std::min(ref.size(), ours.size());
  if (count == 0) {
    err = "No overlapping frames between BLS analyzer and PLUMED file.";
    return false;
  }

  double sumAbsMax = 0.0;
  double sumSqMax = 0.0;
  double sumAbsN = 0.0;
  double sumSqN = 0.0;
  double sumElapsed = 0.0;
  double tauSum = 0.0;
  int tauCount = 0;

  for (std::size_t i = 0; i < count; ++i) {
    const auto& a = ours[i];
    const auto& b = ref[i];
    double diffMax = static_cast<double>(a.maxCluster - b.maxCluster);
    double diffN = static_cast<double>(a.nclusters - b.nclusters);
    double diffElapsed = a.elapsedMs - b.elapsedMs;

    sumAbsMax += std::abs(diffMax);
    sumSqMax += diffMax * diffMax;
    sumAbsN += std::abs(diffN);
    sumSqN += diffN * diffN;
    sumElapsed += diffElapsed;

    if (!a.clusterSizes.empty() && !b.clusterSizes.empty()) {
      tauSum += kendallTau(a.clusterSizes, b.clusterSizes);
      ++tauCount;
    }
  }

  summary.meanAbsMaxCluster = sumAbsMax / static_cast<double>(count);
  summary.rmseMaxCluster = std::sqrt(sumSqMax / static_cast<double>(count));
  summary.meanAbsNClusters = sumAbsN / static_cast<double>(count);
  summary.rmseNClusters = std::sqrt(sumSqN / static_cast<double>(count));
  summary.meanElapsedDiff = sumElapsed / static_cast<double>(count);
  if (tauCount > 0) {
    summary.kendallTau = tauSum / static_cast<double>(tauCount);
  } else {
    summary.kendallTau = 1.0;
  }
  double refTime = 0.0;
  double ourTime = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    refTime += ref[i].elapsedMs;
    ourTime += ours[i].elapsedMs;
  }
  if (ourTime > 0.0) {
    summary.speedup = (refTime > 0.0) ? refTime / ourTime : std::numeric_limits<double>::infinity();
  } else {
    summary.speedup = 0.0;
  }
  summary.matchedFrames = count;
  return true;
}

}  // namespace bls
