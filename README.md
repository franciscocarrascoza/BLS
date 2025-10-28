# BLS Analyzer

Standalone C++ analyzer that implements Bravais Lattice Sampling (BLS) with Skip-DFS
refinement for cluster detection on molecular-dynamics trajectories. The code mirrors
PLUMED structure while staying dependency-light and easy to integrate downstream.

## Features
- Trajectory backends for XTC, TRR, TNG (via optional libraries), GRO, and PDB.
- PLUMED-style configuration file with full control over lattice, connectivity, stride,
  and BLS parameters.
- Per-frame CSV and JSON metrics including cluster counts, seed hits, refined voxels,
  and timing.
- Optional benchmark sidecar with cumulative time and RSS, and comparison utilities
  against PLUMED DFS output.
- Unit tests covering lattice geometry, connectivity modes, PBC handling, and
  determinism.

## Build

```bash
mkdir build && cd build
cmake -DBLS_USE_XDRFILE=ON -DBLS_USE_TNG=ON ..
cmake --build . -j
ctest
```

Disable `BLS_USE_XDRFILE` or `BLS_USE_TNG` if the respective libraries are not
available. OpenMP is detected automatically.

### Docker

```bash
docker build -t bls-analyzer .
docker run --rm bls-analyzer --help
```

## Usage

```bash
./bls_analyze --traj traj.xtc --top topol.gro --conf bls.in \
              --out metrics.csv --json metrics.json --stride 5
```

Optional flags:

- `--bench bench.csv` for cumulative runtime + RSS diagnostics.
- `--compare-plumed colvar.dat` to cross-check against PLUMED DFS output.
- `--threads N` (requires OpenMP build) to enable parallel voxelization.
- `scripts/plot_metrics.sh metrics.csv metrics.png` renders quick-look plots via gnuplot.

### Configuration (`bls.in`)

```
BLS ...
  GROUP ATOMS=all
  BOX AUTO
  GRID_SPACING 0.25
  LATTICE cubic
  CENTERING F
  CONNECTIVITY 6
  DNN 0
  ALPHA 0.7
  RADII 1.5,2.0
  SKIP 3
  OCCUPANCY ANY
  OUTPUT NCLUSTERS,MAX_CLUSTER,SEED_HITS,SEEDS,REFINED_VOXELS
... BLS
```

## Tests

```
ctest --output-on-failure
```

## Output Metrics

Each processed frame adds one row to the CSV report (mirrored in the JSON lines file),
with the columns below:

- `frame`: zero-based frame index.
- `time_ps`: simulation time in picoseconds.
- `natoms`: atoms considered after applying the selection.
- `NX`, `NY`, `NZ`: voxel grid extents along each axis.
- `dNN_vox`: nearest-neighbour spacing in voxel indices used for connectivity.
- `lattice`, `centering`: lattice family and centering active for this run.
- `seeds`: candidate seed voxels produced by the BLS pass.
- `seed_hits`: seeds that yielded occupied voxels after refinement.
- `nclusters`: number of clusters identified in the frame.
- `max_cluster`: size of the largest cluster (voxel count).
- `refined_voxels`: total voxels visited during Skip-DFS refinement.
- `elapsed_ms`: wall-clock milliseconds spent on the frame.

## Status

The analyzer shares core design choices with PLUMED to simplify future integration
while remaining self-contained for standalone benchmarking or deployment.
