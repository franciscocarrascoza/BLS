#include "io/TrajectoryReader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "common/Box.hpp"
#include "util/Logging.hpp"

namespace bls {

class XyzReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    file_.open(path);
    if (!file_) {
      err = "Unable to open XYZ file: " + path;
      return false;
    }
    path_ = path;
    frameIndex_ = 0;
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!file_) {
      err = "XYZ reader not opened.";
      return false;
    }

    std::string line;

    // Line 1: number of atoms
    if (!std::getline(file_, line)) {
      return false;  // EOF, no more frames
    }

    // Skip empty lines between frames
    while (line.empty() || std::all_of(line.begin(), line.end(),
                                       [](unsigned char c) { return std::isspace(c); })) {
      if (!std::getline(file_, line)) {
        return false;
      }
    }

    int natoms = 0;
    try {
      natoms = std::stoi(line);
    } catch (...) {
      err = "Malformed XYZ: invalid atom count on line.";
      return false;
    }

    if (natoms <= 0) {
      err = "Malformed XYZ: atom count must be positive.";
      return false;
    }

    // Line 2: comment line (may contain box info)
    std::string comment;
    if (!std::getline(file_, comment)) {
      err = "Malformed XYZ: missing comment line.";
      return false;
    }

    // Try to parse box info from comment (common formats):
    // Lattice="ax ay az bx by bz cx cy cz" (extended XYZ)
    // or just box dimensions: a b c [alpha beta gamma]
    Mat3 box = parseBoxFromComment(comment);
    double time = parseTimeFromComment(comment);

    // Read atom coordinates
    frame.xyz.clear();
    frame.xyz.reserve(static_cast<std::size_t>(natoms));

    for (int i = 0; i < natoms; ++i) {
      if (!std::getline(file_, line)) {
        err = "Malformed XYZ: premature EOF reading coordinates.";
        return false;
      }

      std::istringstream iss(line);
      std::string element;
      double x, y, z;

      if (!(iss >> element >> x >> y >> z)) {
        err = "Malformed XYZ: invalid coordinate line at atom " + std::to_string(i + 1);
        return false;
      }

      frame.xyz.emplace_back(x, y, z);
    }

    frame.natoms = natoms;
    frame.box = box;
    frame.time = time;

    ++frameIndex_;
    return true;
  }

  void close() override {
    if (file_.is_open()) {
      file_.close();
    }
  }

 private:
  static Mat3 parseBoxFromComment(const std::string& comment) {
    // Try extended XYZ format: Lattice="ax ay az bx by bz cx cy cz"
    auto pos = comment.find("Lattice=\"");
    if (pos != std::string::npos) {
      pos += 9;  // skip 'Lattice="'
      auto endpos = comment.find('"', pos);
      if (endpos != std::string::npos) {
        std::string lattice = comment.substr(pos, endpos - pos);
        std::istringstream iss(lattice);
        double ax, ay, az, bx, by, bz, cx, cy, cz;
        if (iss >> ax >> ay >> az >> bx >> by >> bz >> cx >> cy >> cz) {
          return Mat3{Vec3{ax, ay, az}, Vec3{bx, by, bz}, Vec3{cx, cy, cz}};
        }
      }
    }

    // Try simple box format: "a b c" or "a b c alpha beta gamma"
    std::istringstream iss(comment);
    std::string token;
    std::vector<double> nums;

    while (iss >> token) {
      try {
        double val = std::stod(token);
        nums.push_back(val);
      } catch (...) {
        // Not a number, skip
      }
    }

    if (nums.size() >= 3) {
      double a = nums[0], b = nums[1], c = nums[2];
      if (nums.size() >= 6) {
        // Triclinic: a b c alpha beta gamma
        return buildTriclinicBox(a, b, c, nums[3], nums[4], nums[5]);
      } else {
        // Orthorhombic
        return Mat3{Vec3{a, 0, 0}, Vec3{0, b, 0}, Vec3{0, 0, c}};
      }
    }

    // No box info found, return zero box
    return Mat3{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};
  }

  static double parseTimeFromComment(const std::string& comment) {
    // Look for patterns like "time=X" or "t=X" or "Time: X"
    std::string lower = comment;
    for (auto& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

    for (const char* pattern : {"time=", "t=", "time:", "t:"}) {
      auto pos = lower.find(pattern);
      if (pos != std::string::npos) {
        pos += std::strlen(pattern);
        while (pos < comment.size() && std::isspace(static_cast<unsigned char>(comment[pos]))) {
          ++pos;
        }
        std::string number;
        while (pos < comment.size() &&
               (std::isdigit(static_cast<unsigned char>(comment[pos])) || comment[pos] == '.' ||
                comment[pos] == '-' || comment[pos] == 'e' || comment[pos] == 'E' ||
                comment[pos] == '+')) {
          number.push_back(comment[pos]);
          ++pos;
        }
        if (!number.empty()) {
          try {
            return std::stod(number);
          } catch (...) {
          }
        }
      }
    }

    return 0.0;
  }

  std::ifstream file_;
  std::string path_;
  std::size_t frameIndex_{0};
};

TrajectoryReaderPtr makeXyzReader() { return TrajectoryReaderPtr(new XyzReader()); }

}  // namespace bls
