#include "lattice/Basis.hpp"

#include <algorithm>
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

}  // namespace

double latticeDmin(const Mat3& basis, const std::vector<Vec3>& offsets) {
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

LatticeDescriptor buildLattice(const LatticeSettings& settings) {
  LatticeDescriptor desc;
  desc.basis = buildUnitBasis(settings);
  desc.offsets = centeringOffsets(settings.centering);
  desc.dmin = latticeDmin(desc.basis, desc.offsets);
  return desc;
}

bool commensurateCubicBasis(const LatticeSettings& settings, const LatticeDescriptor& lattice,
                            double dnnVoxel, int nx, int ny, int nz, Mat3& basis,
                            double& effectiveDnn, std::string& err) {
  if (settings.lattice != LatticeType::Cubic) {
    err = "PBC xyz: the probe lattice must be cubic (P, I or F) to be made commensurate with "
          "the periodic cell (audit decision D10); " +
          latticeToString(settings.lattice) + " lattices are not supported under PBC.";
    return false;
  }
  if (!(dnnVoxel > 0.0) || nx < 1 || ny < 1 || nz < 1) {
    err = "PBC xyz: invalid dNN or grid dimensions for the commensurate probe lattice.";
    return false;
  }
  // Unit cubic basis is the identity, so the conventional edge in voxels is dNN / dmin.
  const double a = dnnVoxel / lattice.dmin;
  const int dims[3] = {nx, ny, nz};
  double edge[3];
  for (int k = 0; k < 3; ++k) {
    const int cells = std::max(1, static_cast<int>(std::ceil(dims[k] / a - 1e-9)));
    edge[k] = static_cast<double>(dims[k]) / cells;
  }
  basis = Mat3{Vec3{edge[0], 0.0, 0.0}, Vec3{0.0, edge[1], 0.0}, Vec3{0.0, 0.0, edge[2]}};
  effectiveDnn = latticeDmin(basis, lattice.offsets);
  return true;
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

