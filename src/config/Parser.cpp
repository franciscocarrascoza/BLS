#include "config/Parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "util/Logging.hpp"

namespace bls {

namespace {

std::string trim(const std::string& s) {
  const auto begin =
      std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
  const auto end =
      std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); })
          .base();
  if (begin >= end) {
    return {};
  }
  return std::string(begin, end);
}

std::string toUpper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return s;
}

std::vector<std::string> split(const std::string& line, char delim) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream iss(line);
  while (std::getline(iss, token, delim)) {
    tokens.push_back(trim(token));
  }
  return tokens;
}

std::vector<double> parseDoubleList(const std::string& item) {
  std::vector<double> values;
  auto tokens = split(item, ',');
  values.reserve(tokens.size());
  for (const auto& t : tokens) {
    if (!t.empty()) {
      values.push_back(std::stod(t));
    }
  }
  return values;
}

IndexRange parseRange(const std::string& spec) {
  auto dashPos = spec.find('-');
  if (dashPos == std::string::npos) {
    std::size_t idx = static_cast<std::size_t>(std::stoul(spec));
    return IndexRange{idx - 1, idx - 1};
  }
  std::size_t begin = static_cast<std::size_t>(std::stoul(spec.substr(0, dashPos)));
  std::size_t end = static_cast<std::size_t>(std::stoul(spec.substr(dashPos + 1)));
  if (end < begin) std::swap(begin, end);
  return IndexRange{begin - 1, end - 1};
}

LatticeType latticeFromString(const std::string& s) {
  std::string u = toUpper(s);
  if (u == "CUBIC") return LatticeType::Cubic;
  if (u == "HEXAGONAL") return LatticeType::Hexagonal;
  if (u == "TRICLINIC") return LatticeType::Triclinic;
  throw std::runtime_error("Unsupported lattice type: " + s);
}

CenteringType centeringFromString(const std::string& s) {
  std::string u = toUpper(s);
  if (u == "P") return CenteringType::P;
  if (u == "F") return CenteringType::F;
  if (u == "I") return CenteringType::I;
  throw std::runtime_error("Unsupported centering: " + s);
}

OccupancyMode occupancyFromString(const std::string& s) {
  std::string u = toUpper(s);
  if (u == "ANY") return OccupancyMode::Any;
  if (u == "ALL") return OccupancyMode::All;
  throw std::runtime_error("Unsupported occupancy mode: " + s);
}

}  // namespace

bool Parser::parseFile(const std::string& path, BLSConfig& config, std::string& err) {
  std::ifstream in(path);
  if (!in) {
    err = "Unable to open config file: " + path;
    return false;
  }

  bool inBlock = false;
  std::string line;
  int lineNo = 0;

  while (std::getline(in, line)) {
    ++lineNo;
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;

    if (!inBlock) {
      if (toUpper(line).rfind("BLS", 0) == 0 && line.find("...") != std::string::npos) {
        inBlock = true;
      }
      continue;
    }

    if (line == "... BLS") {
      inBlock = false;
      break;
    }

    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    std::string rest;
    std::getline(iss, rest);
    rest = trim(rest);

    try {
      std::string upperKeyword = toUpper(keyword);
      if (upperKeyword == "GROUP") {
        auto parts = split(rest, '|');
        for (const auto& part : parts) {
          auto eq = part.find('=');
          if (eq == std::string::npos) continue;
          auto key = toUpper(trim(part.substr(0, eq)));
          auto value = trim(part.substr(eq + 1));
          if (key == "ATOMS") {
            auto lower = toUpper(value);
            if (lower == "ALL") {
              config.group.type = GroupSelectorType::All;
              config.group.ranges.clear();
              config.group.names.clear();
            } else if (lower.rfind("INDEX:", 0) == 0) {
              config.group.type = GroupSelectorType::IndexRange;
              config.group.ranges.clear();
              auto rangeSpec = value.substr(6);
              auto tokens = split(rangeSpec, ',');
              for (const auto& tok : tokens) {
                if (!tok.empty()) config.group.ranges.push_back(parseRange(tok));
              }
            } else if (lower.rfind("NAME:", 0) == 0) {
              config.group.type = GroupSelectorType::Name;
              config.group.names = split(value.substr(5), ',');
            }
          }
        }
      } else if (upperKeyword == "BOX") {
        auto tokens = split(rest, ' ');
        if (!tokens.empty() && toUpper(tokens[0]) == "AUTO") {
          config.boxMode = BoxMode::Auto;
        } else {
          config.boxMode = BoxMode::Manual;
          std::unordered_map<std::string, double*> keyMap = {
              {"XLO", &config.manualBox.xlo}, {"XHI", &config.manualBox.xhi},
              {"YLO", &config.manualBox.ylo}, {"YHI", &config.manualBox.yhi},
              {"ZLO", &config.manualBox.zlo}, {"ZHI", &config.manualBox.zhi}};
          for (std::size_t i = 0; i + 1 < tokens.size(); i += 2) {
            auto key = toUpper(tokens[i]);
            auto it = keyMap.find(key);
            if (it != keyMap.end()) {
              *(it->second) = std::stod(tokens[i + 1]);
            }
          }
        }
      } else if (upperKeyword == "GRID_SPACING") {
        config.gridSpacing = std::stod(rest);
      } else if (upperKeyword == "CONNECTIVITY") {
        config.connectivity = std::stoi(rest);
      } else if (upperKeyword == "SKIP") {
        config.skip = std::stoi(rest);
      } else if (upperKeyword == "ALPHA") {
        config.alpha = std::stod(rest);
      } else if (upperKeyword == "DNN") {
        config.dnn = std::stod(rest);
        config.hasExplicitDnn = config.dnn > 0.0;
      } else if (upperKeyword == "RADII") {
        config.radii = parseDoubleList(rest);
      } else if (upperKeyword == "CUTOFF") {
        config.cutoff = std::stod(rest);
      } else if (upperKeyword == "OCCUPANCY") {
        config.occupancy = occupancyFromString(rest);
      } else if (upperKeyword == "LATTICE") {
        config.lattice.lattice = latticeFromString(rest);
      } else if (upperKeyword == "CENTERING") {
        config.lattice.centering = centeringFromString(rest);
      } else if (upperKeyword == "HEX_C_OVER_A") {
        config.lattice.hexCOverA = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_A") {
        config.lattice.triclinicA = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_B") {
        config.lattice.triclinicB = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_C") {
        config.lattice.triclinicC = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_ALPHA") {
        config.lattice.triclinicAlphaDeg = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_BETA") {
        config.lattice.triclinicBetaDeg = std::stod(rest);
      } else if (upperKeyword == "TRICLINIC_GAMMA") {
        config.lattice.triclinicGammaDeg = std::stod(rest);
      } else if (upperKeyword == "STRIDE") {
        config.stride = std::stoi(rest);
      } else if (upperKeyword == "OUTPUT") {
        config.outputs = split(rest, ',');
      } else {
        Logger::warn("Unrecognized keyword at line ", lineNo, ": ", keyword);
      }
    } catch (const std::exception& ex) {
      err = "Parsing error at line " + std::to_string(lineNo) + ": " + ex.what();
      return false;
    }
  }

  if (inBlock) {
    err = "Missing closing \"... BLS\" in config file.";
    return false;
  }

  return true;
}

}  // namespace bls

