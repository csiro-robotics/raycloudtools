## Ray Cloud Tools
A set of command line tools for processing ray clouds, together with an associated C++ library.

1. `$ mkdir build`

2. `$ cd build`

3. `$ cmake ..`

4. `$ make`

To access the tools from anywhere, place in your ~/bashrc:
  export PATH=$PATH:'source code path'/raycloudtools/bin

## Examples:

### rayconvert forest.laz forest_traj.ply
Convert point cloud and trajectory to a single raycloud file:

### raycreate room 0
Generate a single room with a window and door (or other environments), using random seed 0.
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room1.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room3.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
You can visualise the rays in meshlab with Render | Show Vertex Normals. The ray lengths need to be scaled: Tools | Options | NormalLength roughly 0.025 (smaller for larger clouds)
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room2.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raydecimate room.ply 10 cm
Spatially decimates cloud to one point every cubic 10 cm.
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_decimated.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raytranslate room.ply 3 0 0
Translates the ray cloud 3 metres along the x axis.
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_translate.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### rayrotate room.ply 0 0 30
Rotates the ray cloud 30 degrees around the z axis.
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_rotate.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

## raydenoise room.ply 10 cm
Removes rays with isolated end points more than 10 cm from any other, not including unbounded rays.
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise1.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise2.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raysmooth room.ply
Moves ray end points onto the nearest surface, to smooth the resulting cloud.
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_smooth1.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_smooth2.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raytransients min room.ply 0.5 s
Segments out moving or moved objects during the scan, when re-observed more than 0.5 seconds later or before. Leaving the minimum of geometry when transient.

In this raycloud the table and cupboard appear only after the empty room has been scanned for several seconds, so we can isolate these transient objects.
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients1.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
The fixed (untransient) raycloud:
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients2.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>
The remaining transients are also saved out:
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients3.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raycombine all room.ply room2.ply
Combines room and its transformed version together, keeping ***all*** rays.
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_all.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>

### raycombine min room.ply room2.ply
Combines the two ray clouds keeping only the ***minimum*** of geometry where there is a difference. This is a form of union of the two volumes. 
<img img width="640" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_min.png?at=refs%2Fheads%2Fmaster" style="margin-right: 10px;"/>


