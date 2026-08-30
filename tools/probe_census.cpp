// Probe census: the TOTAL number of lattice sites a given lattice places in a
// given grid, independent of occupancy.
//
// E6 compares centerings at equal PROBE BUDGET, so the budget has to be
// measured rather than taken from the continuum formula N = c_L / dNN^3. That
// formula ignores the boundary: the grid is a finite box and the site count is
// an integer lattice point count inside it, which lands a percent or two off
// the continuum value and differently for each centering. Reporting an
// equalization that was only assumed would reproduce exactly the kind of
// unchecked scaling E6 exists to test.
//
// Takes the grid and lattice directly rather than re-deriving them from a
// trajectory: NX, NY, NZ and dNN_vox are all columns of the run's own CSV, so
// the census is computed against precisely the grid that run used.
//
//   bls_probe_census <nx> <ny> <nz> <dnn_voxels> <lattice> <centering>
//   -> "<total_sites> <sites_per_voxel>"
#include <cstdio>
#include <cstdlib>
#include <string>

#include "bls/Options.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"

using namespace bls;

int main(int argc, char** argv) {
  if (argc != 7) {
    std::fprintf(stderr,
                 "usage: %s <nx> <ny> <nz> <dnn_voxels> "
                 "<cubic|hexagonal|triclinic> <P|I|F>\n", argv[0]);
    return 2;
  }
  const int nx = std::atoi(argv[1]);
  const int ny = std::atoi(argv[2]);
  const int nz = std::atoi(argv[3]);
  const double dnn = std::atof(argv[4]);
  const std::string lat = argv[5], cen = argv[6];

  LatticeSettings ls;
  if      (lat == "cubic")      ls.lattice = LatticeType::Cubic;
  else if (lat == "hexagonal")  ls.lattice = LatticeType::Hexagonal;
  else if (lat == "triclinic")  ls.lattice = LatticeType::Triclinic;
  else { std::fprintf(stderr, "unknown lattice '%s'\n", lat.c_str()); return 2; }
  if      (cen == "P") ls.centering = CenteringType::P;
  else if (cen == "I") ls.centering = CenteringType::I;
  else if (cen == "F") ls.centering = CenteringType::F;
  else { std::fprintf(stderr, "unknown centering '%s'\n", cen.c_str()); return 2; }

  const LatticeDescriptor L = buildLattice(ls);
  // Same scaling Analyzer applies: the basis is scaled so the lattice's minimum
  // interpoint distance equals dNN.
  const Mat3 scaled = L.basis * (dnn / L.dmin);

  // The occupancy-free constructor is the volume sweep: every lattice site in
  // the grid, which is the probe budget.
  Enumerator en(scaled, L.offsets, nx, ny, nz);
  const double vol = static_cast<double>(nx) * ny * nz;
  std::printf("%d %.8g\n", en.count(), en.count() / vol);
  return 0;
}
