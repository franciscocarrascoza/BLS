# BLS Analyzer

Standalone C++ analyzer that implements Bravais Lattice Sampling (BLS) with Skip-DFS
refinement for cluster detection on molecular-dynamics trajectories. The code mirrors
PLUMED structure while staying dependency-light and easy to integrate downstream.

## Features
- **Multi-format trajectory support**: XTC, TRR, TNG (via optional libraries), GRO, PDB, XYZ, and MOL/SDF formats.
- **Multiple clustering algorithms**: Choose from 8 different algorithms including BLS (default), traditional DFS, skip-DFS, DBSCAN, hierarchical, k-means, spectral, and GCBD (union-find).
- **PLUMED-style configuration** file with full control over lattice, connectivity, stride, and BLS parameters.
- **Per-frame CSV and JSON metrics** including cluster counts, seed hits, refined voxels, and timing.
- **Optional benchmark sidecar** with cumulative time and RSS, and comparison utilities against PLUMED DFS output.
- **Unit tests** covering lattice geometry, connectivity modes, PBC handling, and determinism.
- **Pure implementations**: All clustering algorithms kept in simple form for fair benchmarking without optimized libraries.

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

### Supported File Formats

**Trajectory files** (via `--traj`):
- **XTC** (GROMACS compressed trajectory)
- **TRR** (GROMACS full-precision trajectory)
- **TNG** (GROMACS next-generation format, requires TNG library)
- **GRO** (GROMACS structure file)
- **PDB** (Protein Data Bank format, supports multi-MODEL trajectories)
- **XYZ** (Extended XYZ format with optional lattice information)
- **MOL/SDF** (MDL Molfile format, V2000 and V3000)

Use `--format xyz` or `--format mol` to override automatic format detection.

### Clustering Algorithms

Select a clustering algorithm using `--algo ALGORITHM`:

- **`bls`** (default) - Bravais Lattice Sampling with Skip-DFS refinement
- **`traditional_dfs`** - Traditional depth-first search clustering
- **`skip_dfs`** - Skip-DFS without lattice sampling
- **`dbscan`** - DBSCAN with grid-based spatial indexing
- **`hierarchical`** - Single-linkage hierarchical clustering
- **`kmeans`** - K-means clustering
- **`spectral`** - Simplified spectral clustering
- **`gcbd`** - Grid-based connectivity using union-find

**Algorithm-specific parameters:**
- `--algo-skip N` - Skip distance for skip_dfs (default: 3)
- `--algo-eps F` - Epsilon for DBSCAN (default: 3.0)
- `--algo-minpts N` - MinPts for DBSCAN (default: 10)
- `--algo-k N` - K for k-means/spectral (default: 20)
- `--algo-threshold F` - Threshold for hierarchical (default: 4.0)

**Example with algorithm selection:**
```bash
# Use GCBD (fast union-find clustering)
./bls_analyze --traj validation/g1000_20Ih.pdb --conf bls.in --algo gcbd

# Use DBSCAN with custom parameters
./bls_analyze --traj traj.xyz --conf bls.in --algo dbscan --algo-eps 2.5 --algo-minpts 15

# Use traditional DFS
./bls_analyze --traj traj.xtc --conf bls.in --algo traditional_dfs
```

### Optional Flags

- `--bench bench.csv` for cumulative runtime + RSS diagnostics.
- `--compare-plumed colvar.dat` to cross-check against PLUMED DFS output.
- `--threads N` (requires OpenMP build) to enable parallel voxelization.
- `scripts/plot_metrics.sh metrics.csv metrics.png` renders quick-look plots via gnuplot.

### Configuration File (`bls.in`)

The configuration file uses a PLUMED-style format. All parameters are optional and have sensible defaults:

```
BLS ...
  # Atom selection
  GROUP ATOMS=all                    # 'all' or index ranges like '1-100' or names like 'name:CA'

  # Box configuration
  BOX AUTO                           # AUTO (use trajectory box) or MANUAL with xlo,xhi,ylo,yhi,zlo,zhi

  # Grid parameters
  GRID_SPACING 0.25                  # Voxel size in Angstroms
  CONNECTIVITY 6                     # 6 (face neighbors) or 26 (face+edge+corner)

  # BLS lattice parameters (only used when --algo bls)
  LATTICE cubic                      # cubic, hexagonal, or triclinic
  CENTERING F                        # P (primitive), F (face-centered), or I (body-centered)
  DNN 0                              # Nearest-neighbor distance (0 = auto-compute)
  ALPHA 0.7                          # Scaling factor for auto dNN
  RADII 1.5,2.0                      # Atomic radii for auto dNN computation

  # Refinement parameters
  SKIP 3                             # Skip distance for Skip-DFS

  # Rasterization
  CUTOFF 0.6                         # Cutoff radius for voxel occupancy (Angstroms)
  OCCUPANCY ANY                      # ANY (atom touches voxel) or ALL (atom fully inside)

  # Output selection
  OUTPUT NCLUSTERS,MAX_CLUSTER,SEED_HITS,SEEDS,REFINED_VOXELS
... BLS
```

A fully commented example is available in `validation/bls_test.in`.

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
- `dNN_vox`: nearest-neighbour spacing in voxel indices used for connectivity (BLS only; 0 for other algorithms).
- `lattice`, `centering`: lattice family and centering active for this run (for BLS) or algorithm name (for other algorithms).
- `seeds`: candidate seed voxels produced by the BLS pass (0 for non-BLS algorithms).
- `seed_hits`: seeds that yielded occupied voxels after refinement (0 for non-BLS algorithms).
- `nclusters`: number of clusters identified in the frame.
- `max_cluster`: size of the largest cluster (voxel count).
- `refined_voxels`: total voxels visited during clustering.
- `elapsed_ms`: wall-clock milliseconds spent on the frame.

## Algorithm Benchmarking

All clustering algorithms are implemented in their pure, unoptimized form to enable fair performance comparisons for research purposes. This means:

- **No external optimization libraries** - All algorithms use only standard C++ without BLAS, Eigen, or similar optimized libraries
- **Simple, readable implementations** - Code follows the textbook algorithms directly for proof-of-concept evaluation
- **Independent modules** - Each algorithm is self-contained in `src/cluster/Algorithms.cpp` for easy inspection

This design allows researchers to:
1. Compare time complexity of different clustering approaches on the same grid data
2. Evaluate trade-offs between algorithm sophistication and runtime
3. Benchmark BLS lattice sampling against traditional clustering methods
4. Maintain reproducible, transparent algorithm implementations

For production use cases requiring maximum performance, consider:
- Using optimized library implementations (scipy, sklearn, etc.) in Python
- Enabling compiler optimizations (`-O3 -march=native`)
- Profiling and optimizing the hot paths specific to your workload

## Status

The analyzer shares core design choices with PLUMED to simplify future integration
while remaining self-contained for standalone benchmarking or deployment.
