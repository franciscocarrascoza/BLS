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
        // BOX AUTO | BOX CELL | BOX [MANUAL] XLO v XHI v YLO v YHI v ZLO v ZHI v
        // The manual form used to accept anything: a leading MANUAL token shifted the
        // key/value pairing so every bound silently stayed 0, and unknown keys were
        // ignored (audit HARN-2). All six bounds are now required and checked.
        std::vector<std::string> tokens;
        {
          std::istringstream ts(rest);
          std::string t;
          while (ts >> t) tokens.push_back(t);
        }
        const std::string head = tokens.empty() ? std::string() : toUpper(tokens[0]);
        if (head == "AUTO" && tokens.size() == 1) {
          config.boxMode = BoxMode::Auto;
        } else if (head == "CELL" && tokens.size() == 1) {
          config.boxMode = BoxMode::Cell;
        } else {
          config.boxMode = BoxMode::Manual;
          std::unordered_map<std::string, double*> keyMap = {
              {"XLO", &config.manualBox.xlo}, {"XHI", &config.manualBox.xhi},
              {"YLO", &config.manualBox.ylo}, {"YHI", &config.manualBox.yhi},
              {"ZLO", &config.manualBox.zlo}, {"ZHI", &config.manualBox.zhi}};
          std::size_t first = (head == "MANUAL") ? 1 : 0;
          if ((tokens.size() - first) != 12) {
            throw std::runtime_error(
                "BOX expects AUTO, CELL, or [MANUAL] XLO v XHI v YLO v YHI v ZLO v ZHI v");
          }
          std::unordered_map<std::string, bool> seen;
          for (std::size_t i = first; i + 1 < tokens.size(); i += 2) {
            auto key = toUpper(tokens[i]);
            auto it = keyMap.find(key);
            if (it == keyMap.end() || seen[key]) {
              throw std::runtime_error("BOX: unknown or repeated bound '" + tokens[i] + "'");
            }
            seen[key] = true;
            *(it->second) = std::stod(tokens[i + 1]);
          }
        }
      } else if (upperKeyword == "PBC") {
        const std::string v = toUpper(trim(rest));
        if (v == "NONE") {
          config.pbc = PeriodicAxes{};
        } else if (v == "XYZ") {
          config.pbc = PeriodicAxes{true, true, true};
        } else {
          throw std::runtime_error("PBC expects 'xyz' or 'none' (per-axis periodicity is not "
                                   "implemented by any method)");
        }
      } else if (upperKeyword == "LATTICE_ORIGIN") {
        std::istringstream ts(rest);
        double o[3];
        if (!(ts >> o[0] >> o[1] >> o[2])) {
          throw std::runtime_error("LATTICE_ORIGIN expects three numbers (voxels)");
        }
        for (int k = 0; k < 3; ++k) config.latticeOrigin[k] = o[k];
      } else if (upperKeyword == "GRID_SPACING") {
        config.gridSpacing = std::stod(rest);
      } else if (upperKeyword == "CONNECTIVITY") {
        config.connectivity = std::stoi(rest);
      } else if (upperKeyword == "SKIP") {
        // Kept as the config spelling for backward compatibility with the 83
        // E0-E5 decks. It maps to the refinement stride and to nothing else --
        // not to --algo-skip, which belongs to cluster::skipDFS.
        config.refinementStride = std::stoi(rest);
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
        config.lattice.latticeSet = true;
      } else if (upperKeyword == "CENTERING") {
        config.lattice.centering = centeringFromString(rest);
        config.lattice.centeringSet = true;
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

  // Deliberately fatal rather than defaulted. See LatticeSettings in
  // Options.hpp: the old silent cubic/F fallback shifts the cluster count by
  // 17% on ld-asw now that lattice enumeration actually works, and a config
  // that does not say which lattice it means cannot be reproduced from its
  // own text.
  if (!config.lattice.latticeSet || !config.lattice.centeringSet) {
    std::string missing;
    if (!config.lattice.latticeSet) missing = "LATTICE";
    if (!config.lattice.centeringSet) {
      if (!missing.empty()) missing += " and ";
      missing += "CENTERING";
    }
    err = "Config is missing " + missing +
          ". Both are required -- there is no default. State the lattice the "
          "run is meant to use, e.g.\n"
          "    LATTICE cubic\n"
          "    CENTERING F\n"
          "(LATTICE: cubic | hexagonal | triclinic. CENTERING: P | F | I.)";
    return false;
  }

  if (config.pbc.any() && config.boxMode == BoxMode::Auto) {
    err = "PBC xyz needs a periodic cell: use BOX CELL (the file's CRYST1 cell) or BOX "
          "MANUAL. BOX AUTO may synthesise a bounding box around the atoms, which is not a "
          "periodic cell.";
    return false;
  }

  return true;
}

}  // namespace bls

