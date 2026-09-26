<!-- Created: 2026-09-26T18:15+02:00 | Rewritten: 2026-09-26T18:34+02:00 | by: CC audit (manager) | purpose: sources of third-party methods and parameter values used by the comparison algorithms (finding VCCS-4) -->
# Third-party notices and sources

## VCCS (`vccs_optimized`, `src/cluster/Algorithms.cpp`)

Implemented from the publication, not from any library's code: J. Papon, A. Abramov,
M. Schoeler, F. Wörgötter, "Voxel Cloud Connectivity Segmentation — Supervoxels for Point
Clouds", IEEE CVPR 2013, doi:10.1109/CVPR.2013.264 (§3.1–3.4; the code
comments quote the passages each rule follows).

Two values the paper does not state are taken from the Point Cloud Library's
`SupervoxelClustering` (`segmentation/include/pcl/segmentation/impl/supervoxel_clustering.hpp`,
written by the paper's first author; PCL is distributed under the BSD licence,
www.pointclouds.org): the seed-filter search radius R_search = 0.5·R_seed, and the count of
voxels "occupied by a planar surface" through the search volume taken as its area over the
voxel area, πR_search²/R_voxel². PCL multiplies that count by 0.05, noting it is "smaller than
the value mentioned in the original paper"; that factor is not used here. No PCL source code
is included in this repository.
