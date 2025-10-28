#include "io/TrajectoryReader.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "util/Logging.hpp"

namespace bls {

class GroReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    file_.open(path);
    if (!file_) {
      err = "Unable to open GRO file: " + path;
      return false;
    }
    path_ = path;
    frameIndex_ = 0;
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!file_) {
      err = "GRO reader not opened.";
      return false;
    }

    std::string title;
    if (!std::getline(file_, title)) {
      return false;
    }

    std::string natomLine;
    if (!std::getline(file_, natomLine)) {
      err = "Malformed GRO: missing atom count line.";
      return false;
    }
    int natoms = std::stoi(natomLine);
    frame.xyz.clear();
    frame.xyz.reserve(static_cast<std::size_t>(natoms));

    std::string line;
    for (int i = 0; i < natoms; ++i) {
      if (!std::getline(file_, line)) {
        err = "Malformed GRO: premature EOF for coordinates.";
        return false;
      }
      if (line.size() < 44) {
        err = "Malformed GRO: coordinate line too short.";
        return false;
      }
      double x = std::stod(line.substr(20, 8));
      double y = std::stod(line.substr(28, 8));
      double z = std::stod(line.substr(36, 8));
      frame.xyz.emplace_back(x, y, z);
    }

    if (!std::getline(file_, line)) {
      err = "Malformed GRO: missing box line.";
      return false;
    }

    std::istringstream boxss(line);
    double x = 0, y = 0, z = 0;
    boxss >> x >> y >> z;
    frame.box = Mat3{Vec3{x, 0, 0}, Vec3{0, y, 0}, Vec3{0, 0, z}};
    frame.natoms = natoms;
    frame.time = extractTime(title);

    ++frameIndex_;
    return true;
  }

  void close() override {
    if (file_.is_open()) {
      file_.close();
    }
  }

 private:
  static double extractTime(const std::string& header) {
    auto pos = header.find("t=");
    if (pos == std::string::npos) return 0.0;
    pos += 2;
    std::string number;
    while (pos < header.size() && std::isspace(static_cast<unsigned char>(header[pos]))) ++pos;
    while (pos < header.size() &&
           (std::isdigit(static_cast<unsigned char>(header[pos])) || header[pos] == '.' ||
            header[pos] == '-')) {
      number.push_back(header[pos]);
      ++pos;
    }
    if (number.empty()) return 0.0;
    return std::stod(number);
  }

  std::ifstream file_;
  std::string path_;
  std::size_t frameIndex_{0};
};

TrajectoryReaderPtr makeGroReader() { return TrajectoryReaderPtr(new GroReader()); }

}  // namespace bls

