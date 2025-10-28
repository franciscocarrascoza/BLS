#include "io/Topology.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace bls {

namespace {

std::string trim(const std::string& s) {
  const auto begin =
      std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
  const auto end =
      std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); })
          .base();
  if (begin >= end) return {};
  return std::string(begin, end);
}

bool loadGro(const std::string& path, Topology& topo, std::string& err) {
  std::ifstream in(path);
  if (!in) {
    err = "Unable to open GRO topology: " + path;
    return false;
  }
  std::string line;
  if (!std::getline(in, line)) {
    err = "Empty GRO file.";
    return false;
  }
  if (!std::getline(in, line)) {
    err = "Malformed GRO: missing atom count.";
    return false;
  }
  int natoms = std::stoi(line);
  topo.atoms.clear();
  topo.atoms.reserve(static_cast<std::size_t>(natoms));
  for (int i = 0; i < natoms; ++i) {
    if (!std::getline(in, line)) {
      err = "Malformed GRO: unexpected EOF while reading atoms.";
      return false;
    }
    if (line.size() < 15) continue;
    std::string name = trim(line.substr(10, 5));
    topo.atoms.push_back(AtomInfo{name, i});
  }
  return true;
}

bool loadPdb(const std::string& path, Topology& topo, std::string& err) {
  std::ifstream in(path);
  if (!in) {
    err = "Unable to open PDB topology: " + path;
    return false;
  }
  topo.atoms.clear();
  std::string line;
  int index = 0;
  while (std::getline(in, line)) {
    if (line.size() < 6) continue;
    std::string rec = line.substr(0, 6);
    std::transform(rec.begin(), rec.end(), rec.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (rec == "ATOM  " || rec == "HETATM") {
      std::string name = line.size() >= 16 ? trim(line.substr(12, 4)) : "";
      topo.atoms.push_back(AtomInfo{name, index++});
    }
  }
  if (topo.atoms.empty()) {
    err = "No atoms found in PDB topology.";
    return false;
  }
  return true;
}

}  // namespace

bool loadTopology(const std::string& path, Topology& topo, std::string& err) {
  auto pos = path.find_last_of('.');
  std::string ext = (pos == std::string::npos) ? "" : path.substr(pos + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (ext == "gro") {
    return loadGro(path, topo, err);
  }
  if (ext == "pdb") {
    return loadPdb(path, topo, err);
  }
  err = "Unsupported topology format: " + ext;
  return false;
}

}  // namespace bls

