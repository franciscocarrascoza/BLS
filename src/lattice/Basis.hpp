#pragma once

#include <string>
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

// Smallest distance between two sites of the lattice basis*(idx + offset).
double latticeDmin(const Mat3& basis, const std::vector<Vec3>& offsets);

// Periodic grids (deck PBC xyz; audit decision D10). The probe lattice must tile the
// periodic N_x x N_y x N_z voxel grid, or sites would crowd or thin out across the
// boundary. Lengths are in voxels of the grid the lattice is rounded onto. The
// conventional cubic cell edge a = dnnVoxel / dmin (unit basis) is replaced per axis by
//   n_i = ceil(N_i / a) cells,  a_i = N_i / n_i <= a,
// so the lattice is never sparser than requested. `basis` = diag(a_x, a_y, a_z);
// `effectiveDnn` = its nearest-neighbour distance in voxels (<= dnnVoxel), the value to
// report. Only cubic lattices (P, I, F) have this form; hexagonal and triclinic probe
// lattices under PBC are rejected (err set, false returned).
bool commensurateCubicBasis(const LatticeSettings& settings, const LatticeDescriptor& lattice,
                            double dnnVoxel, int nx, int ny, int nz, Mat3& basis,
                            double& effectiveDnn, std::string& err);

std::string latticeToString(LatticeType lattice);
std::string centeringToString(CenteringType centering);

}  // namespace bls

