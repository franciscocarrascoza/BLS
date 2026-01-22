#include "io/TrajectoryReader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "util/Logging.hpp"

namespace bls {

class MolReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    file_.open(path);
    if (!file_) {
      err = "Unable to open MOL/SDF file: " + path;
      return false;
    }
    path_ = path;
    frameIndex_ = 0;
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!file_) {
      err = "MOL reader not opened.";
      return false;
    }

    std::string line;

    // Skip to next molecule (for SDF files with multiple records)
    // or start reading from beginning

    // Header block: 3 lines (name, program/timestamp, comment)
    std::string molName;
    if (!std::getline(file_, molName)) {
      return false;  // EOF, no more frames
    }

    // Skip empty lines or "$$$$" markers from previous molecule
    while (molName.empty() || molName.find("$$$$") != std::string::npos ||
           std::all_of(molName.begin(), molName.end(),
                       [](unsigned char c) { return std::isspace(c); })) {
      if (!std::getline(file_, molName)) {
        return false;
      }
    }

    // Check if this is a V3000 format (starts with specific marker)
    bool isV3000 = (molName.find("V3000") != std::string::npos);

    // Line 2: program/timestamp info
    std::string programLine;
    if (!std::getline(file_, programLine)) {
      err = "Malformed MOL: missing program line.";
      return false;
    }

    // Line 3: comment
    std::string commentLine;
    if (!std::getline(file_, commentLine)) {
      err = "Malformed MOL: missing comment line.";
      return false;
    }

    // Line 4: counts line
    std::string countsLine;
    if (!std::getline(file_, countsLine)) {
      err = "Malformed MOL: missing counts line.";
      return false;
    }

    // Check for V3000 block indicator
    if (countsLine.find("V3000") != std::string::npos) {
      return readV3000(frame, err);
    }

    // V2000 format: counts line "aaabbblllfff..."
    // First 3 chars: atom count, next 3: bond count
    int natoms = 0;
    int nbonds = 0;

    try {
      if (countsLine.size() >= 3) {
        natoms = std::stoi(countsLine.substr(0, 3));
      }
      if (countsLine.size() >= 6) {
        nbonds = std::stoi(countsLine.substr(3, 3));
      }
    } catch (...) {
      // Try whitespace-separated parsing
      std::istringstream iss(countsLine);
      if (!(iss >> natoms >> nbonds)) {
        err = "Malformed MOL: invalid counts line.";
        return false;
      }
    }

    if (natoms <= 0) {
      err = "Malformed MOL: atom count must be positive.";
      return false;
    }

    // Read atom block
    frame.xyz.clear();
    frame.xyz.reserve(static_cast<std::size_t>(natoms));

    for (int i = 0; i < natoms; ++i) {
      if (!std::getline(file_, line)) {
        err = "Malformed MOL: premature EOF reading atom block.";
        return false;
      }

      double x, y, z;
      // V2000 format: columns 0-10 (x), 10-20 (y), 20-30 (z)
      // Each coordinate is 10 characters wide
      if (line.size() >= 30) {
        try {
          x = std::stod(line.substr(0, 10));
          y = std::stod(line.substr(10, 10));
          z = std::stod(line.substr(20, 10));
        } catch (...) {
          // Try whitespace-separated
          std::istringstream iss(line);
          if (!(iss >> x >> y >> z)) {
            err = "Malformed MOL: invalid coordinate at atom " + std::to_string(i + 1);
            return false;
          }
        }
      } else {
        std::istringstream iss(line);
        if (!(iss >> x >> y >> z)) {
          err = "Malformed MOL: invalid coordinate at atom " + std::to_string(i + 1);
          return false;
        }
      }

      frame.xyz.emplace_back(x, y, z);
    }

    // Skip bond block and properties until end of record
    while (std::getline(file_, line)) {
      if (line.find("M  END") != std::string::npos) {
        break;
      }
      if (line.find("$$$$") != std::string::npos) {
        break;
      }
    }

    frame.natoms = natoms;
    frame.box = Mat3{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};  // MOL files don't have box info
    frame.time = static_cast<double>(frameIndex_);

    ++frameIndex_;
    return true;
  }

  void close() override {
    if (file_.is_open()) {
      file_.close();
    }
  }

 private:
  bool readV3000(Frame& frame, std::string& err) {
    std::string line;
    int natoms = 0;

    // Find BEGIN CTAB
    while (std::getline(file_, line)) {
      if (line.find("BEGIN CTAB") != std::string::npos) {
        break;
      }
    }

    // Find COUNTS line
    while (std::getline(file_, line)) {
      if (line.find("COUNTS") != std::string::npos) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
          if (token == "COUNTS") continue;
          if (token == "M" || token == "V30") continue;
          try {
            natoms = std::stoi(token);
            break;
          } catch (...) {
          }
        }
        break;
      }
    }

    if (natoms <= 0) {
      err = "Malformed V3000 MOL: could not parse atom count.";
      return false;
    }

    // Find BEGIN ATOM
    while (std::getline(file_, line)) {
      if (line.find("BEGIN ATOM") != std::string::npos) {
        break;
      }
    }

    // Read atoms
    frame.xyz.clear();
    frame.xyz.reserve(static_cast<std::size_t>(natoms));

    for (int i = 0; i < natoms; ++i) {
      if (!std::getline(file_, line)) {
        err = "Malformed V3000 MOL: premature EOF reading atoms.";
        return false;
      }

      if (line.find("END ATOM") != std::string::npos) {
        break;
      }

      // V3000 atom line: M  V30 index type x y z ...
      std::istringstream iss(line);
      std::string token;
      int fieldIdx = 0;
      double x = 0, y = 0, z = 0;

      while (iss >> token) {
        if (token == "M" || token == "V30") continue;
        fieldIdx++;
        if (fieldIdx == 1) continue;       // index
        if (fieldIdx == 2) continue;       // atom type
        if (fieldIdx == 3) x = std::stod(token);
        if (fieldIdx == 4) y = std::stod(token);
        if (fieldIdx == 5) {
          z = std::stod(token);
          break;
        }
      }

      frame.xyz.emplace_back(x, y, z);
    }

    // Skip to end of record
    while (std::getline(file_, line)) {
      if (line.find("M  END") != std::string::npos) {
        break;
      }
      if (line.find("$$$$") != std::string::npos) {
        break;
      }
    }

    frame.natoms = static_cast<int>(frame.xyz.size());
    frame.box = Mat3{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};
    frame.time = static_cast<double>(frameIndex_);

    ++frameIndex_;
    return true;
  }

  std::ifstream file_;
  std::string path_;
  std::size_t frameIndex_{0};
};

TrajectoryReaderPtr makeMolReader() { return TrajectoryReaderPtr(new MolReader()); }

}  // namespace bls
