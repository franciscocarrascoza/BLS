<!-- Created: 2026-09-26T18:15+02:00 | by: CC audit (manager) | purpose: third-party notices for code derived from other projects (brief: PCL-derived code keeps its BSD-3 notice; finding VCCS-4) -->
# Third-party notices

## Point Cloud Library (PCL) — supervoxel clustering (VCCS)

`src/cluster/Algorithms.cpp` (`vccsOptimized`) ports the supervoxel algorithm of PCL's
`segmentation/include/pcl/segmentation/impl/supervoxel_clustering.hpp` (seed selection and
pruning, `expandSupervoxels`, `SupervoxelHelper::expand`, `updateCentroid`) onto BLS's occupancy
grid, keeping only the spatial distance term. Method: J. Papon, A. Abramov, M. Schoeler,
F. Wörgötter, "Voxel Cloud Connectivity Segmentation — Supervoxels for Point Clouds", CVPR 2013,
doi:10.1109/CVPR.2013.264. The PCL notice, reproduced from the header of that file:

```
Software License Agreement (BSD License)

 Point Cloud Library (PCL) - www.pointclouds.org

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.
  * Neither the name of Willow Garage, Inc. nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.

Author : jpapon@gmail.com
```

The neighbour-count rule of the seed pruning follows FLANN's radius search (strict `dist < radius`),
which PCL's kd-tree search uses (flann-lib/flann, BSD licence); no FLANN code is included.
