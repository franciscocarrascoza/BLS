#include "lattice/Basis.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/Box.hpp"
#include "common/Types.hpp"

namespace bls {

namespace {

Mat3 buildUnitBasis(const LatticeSettings& settings) {
  switch (settings.lattice) {
    case LatticeType::Cubic: {
      return Mat3{Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    }
    case LatticeType::Hexagonal: {
      const double sqrt3 = std::sqrt(3.0);
      Vec3 a1{1.0, 0.0, 0.0};
      Vec3 a2{-0.5, 0.5 * sqrt3, 0.0};
      Vec3 a3{0.0, 0.0, settings.hexCOverA};
      return Mat3{a1, a2, a3};
    }
    case LatticeType::Triclinic: {
      return buildTriclinicBox(settings.triclinicA, settings.triclinicB, settings.triclinicC,
                               settings.triclinicAlphaDeg, settings.triclinicBetaDeg,
                               settings.triclinicGammaDeg);
    }
  }
  throw std::runtime_error("Unhandled lattice type");
}

std::vector<Vec3> centeringOffsets(CenteringType centering) {
  if (centering == CenteringType::P) {
    return {Vec3{0.0, 0.0, 0.0}};
  }
  if (centering == CenteringType::F) {
    return {Vec3{0.0, 0.0, 0.0}, Vec3{0.5, 0.5, 0.0}, Vec3{0.5, 0.0, 0.5}, Vec3{0.0, 0.5, 0.5}};
  }
  if (centering == CenteringType::I) {
    return {Vec3{0.0, 0.0, 0.0}, Vec3{0.5, 0.5, 0.5}};
  }
  throw std::runtime_error("Unhandled centering type");
}

double computeDmin(const Mat3& basis, const std::vector<Vec3>& offsets) {
  double dmin = std::numeric_limits<double>::infinity();
  const int range = 1;
  for (std::size_t i = 0; i < offsets.size(); ++i) {
    Vec3 vi = offsets[i];
    for (int ix = -range; ix <= range; ++ix) {
      for (int iy = -range; iy <= range; ++iy) {
        for (int iz = -range; iz <= range; ++iz) {
          Vec3 shift{static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)};
          for (std::size_t j = 0; j < offsets.size(); ++j) {
            if (ix == 0 && iy == 0 && iz == 0 && i == j) continue;
            Vec3 vj = offsets[j] + shift;
            Vec3 diff = basis * (vj - vi);
            double dist = norm(diff);
            if (dist > 1e-8 && dist < dmin) {
              dmin = dist;
            }
          }
        }
      }
    }
  }
  return dmin;
}

}  // namespace

LatticeDescriptor buildLattice(const LatticeSettings& settings) {
  LatticeDescriptor desc;
  desc.basis = buildUnitBasis(settings);
  desc.offsets = centeringOffsets(settings.centering);
  desc.dmin = computeDmin(desc.basis, desc.offsets);
  return desc;
}

std::string latticeToString(LatticeType lattice) {
  switch (lattice) {
    case LatticeType::Cubic:
      return "cubic";
    case LatticeType::Hexagonal:
      return "hexagonal";
    case LatticeType::Triclinic:
      return "triclinic";
  }
  return "unknown";
}

std::string centeringToString(CenteringType centering) {
  switch (centering) {
    case CenteringType::P:
      return "P";
    case CenteringType::F:
      return "F";
    case CenteringType::I:
      return "I";
  }
  return "?";
}

}  // namespace bls

