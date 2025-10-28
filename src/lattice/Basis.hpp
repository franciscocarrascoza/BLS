#pragma once

#include <vector>

#include "bls/Options.hpp"
#include "common/Types.hpp"

namespace bls {

struct LatticeDescriptor {
  Mat3 basis;
  std::vector<Vec3> offsets;
  double dmin;
};

LatticeDescriptor buildLattice(const LatticeSettings& settings);

std::string latticeToString(LatticeType lattice);
std::string centeringToString(CenteringType centering);

}  // namespace bls

