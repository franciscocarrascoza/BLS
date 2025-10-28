#pragma once

#include <string>
#include <vector>

namespace bls {

struct AtomInfo {
  std::string name;
  int index{0};
};

struct Topology {
  std::vector<AtomInfo> atoms;
};

bool loadTopology(const std::string& path, Topology& topo, std::string& err);

}  // namespace bls

