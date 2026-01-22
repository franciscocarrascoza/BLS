#include "io/TrajectoryReader.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "common/Box.hpp"
#include "util/Logging.hpp"

namespace bls {

class PdbReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    file_.open(path);
    if (!file_) {
      err = "Unable to open PDB file: " + path;
      return false;
    }
    path_ = path;
    parsedCrystal_ = false;
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!file_.is_open()) {
      err = "PDB reader not opened.";
      return false;
    }

    // Check for EOF from previous read
    if (file_.eof()) {
      return false;  // Normal EOF, no error
    }

    std::string line;
    std::vector<Vec3> coords;
    double modelTime = 0.0;
    bool inModel = false;

    while (std::getline(file_, line)) {
      if (line.size() < 6) continue;
      std::string rec = line.substr(0, 6);
      for (auto& c : rec) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

      if (rec == "CRYST1") {
        if (!parsedCrystal_) {
          parseCrystal(line);
        }
      } else if (rec == "MODEL ") {
        if (inModel && !coords.empty()) {
          break;
        }
        inModel = true;
        modelTime = parseModelTime(line);
        coords.clear();
      } else if (rec == "ENDMDL") {
        if (!coords.empty()) {
          break;
        }
      } else if (rec == "ATOM  " || rec == "HETATM") {
        if (line.size() < 54) continue;
        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));
        coords.emplace_back(x, y, z);
      }
    }

    if (coords.empty()) {
      return false;
    }

    frame.xyz = std::move(coords);
    frame.natoms = static_cast<int>(frame.xyz.size());
    frame.time = modelTime;
    frame.box = box_;

    return true;
  }

  void close() override {
    if (file_.is_open()) {
      file_.close();
    }
  }

 private:
  void parseCrystal(const std::string& line) {
    if (line.size() < 66) return;
    double a = std::stod(line.substr(6, 9));
    double b = std::stod(line.substr(15, 9));
    double c = std::stod(line.substr(24, 9));
    double alpha = std::stod(line.substr(33, 7));
    double beta = std::stod(line.substr(40, 7));
    double gamma = std::stod(line.substr(47, 7));
    box_ = buildTriclinicBox(a, b, c, alpha, beta, gamma);
    parsedCrystal_ = true;
  }

  static double parseModelTime(const std::string& line) {
    auto pos = line.find("TIME=");
    if (pos == std::string::npos) return 0.0;
    pos += 5;
    std::string number;
    while (pos < line.size() &&
           (std::isspace(static_cast<unsigned char>(line[pos])) || line[pos] == '=')) {
      ++pos;
    }
    while (pos < line.size() &&
           (std::isdigit(static_cast<unsigned char>(line[pos])) || line[pos] == '.' ||
            line[pos] == '-')) {
      number.push_back(line[pos]);
      ++pos;
    }
    if (number.empty()) return 0.0;
    return std::stod(number);
  }

  std::ifstream file_;
  std::string path_;
  Mat3 box_{Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}};
  bool parsedCrystal_{false};
};

TrajectoryReaderPtr makePdbReader() { return TrajectoryReaderPtr(new PdbReader()); }

}  // namespace bls

