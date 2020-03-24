## Ray Cloud Tools
A set of command line tools for processing ray clouds, together with an associated C++ library.

1. `$ mkdir build`

2. `$ cd build`

3. `$ cmake ..`

4. `$ make`

To access the tools from anywhere, place in your ~/bashrc:
  export PATH=$PATH:'source code path'/raycloudtools/bin

## Examples:

**rayconvert forest.laz forest_traj.ply**  Convert point cloud and trajectory to a single raycloud file: forest.ply

**raycreate room 0**  Generate a single room with a window and door (or other environments), using random seed 0.
<p align="center">
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room1.png?at=refs%2Fheads%2Fmaster"/>
  <img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room3.png?at=refs%2Fheads%2Fmaster"/>
</p>
You can visualise the rays in meshlab with Render | Show Vertex Normals. The ray lengths need to be scaled: Tools | Options | NormalLength roughly 0.025 (smaller for larger clouds)

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room2.png?at=refs%2Fheads%2Fmaster"/></p>

**raydecimate room.ply 10 cm**  Spatially decimate cloud to one point every cubic 10 cm.

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_decimated.png?at=refs%2Fheads%2Fmaster"/></p>

**raytranslate room.ply 3 0 0**  Translate the ray cloud 3 metres along the x axis.

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_translate.png?at=refs%2Fheads%2Fmaster"/></p>

**rayrotate room.ply 0 0 30**  Rotate the ray cloud 30 degrees around the z axis.

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_rotate.png?at=refs%2Fheads%2Fmaster"/></p>

**raydenoise room.ply 10 cm**  Remove rays with isolated end points more than 10 cm from any other, not including unbounded rays.

<p align="center">
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise2.png?at=refs%2Fheads%2Fmaster"/>

**raysmooth room.ply**  Move ray end points onto the nearest surface, to smooth the resulting cloud.

<p align="center">
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_smooth1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_smooth2.png?at=refs%2Fheads%2Fmaster"/>
</p>

**raytransients min room.ply 0.5 s**  Segment out moving or moved objects during the scan, when re-observed more than 0.5 seconds later or before. 

Leaving the minimum of geometry when transient.

In this raycloud the table and cupboard appear only after the empty room has been scanned for several seconds, so we can isolate these transient objects.

<p align="center">
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients2.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_transients3.png?at=refs%2Fheads%2Fmaster"/>
</p>

         Original cloud                                    The fixed (untransient) raycloud                             the remaining transient rays are also saved out

**raycombine all room.ply room2.ply**  Combine room and its transformed version together, keeping ***all*** rays.

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_all.png?at=refs%2Fheads%2Fmaster"/></p>

**raycombine min room.ply room2.ply**  Combine the two ray clouds keeping only the ***minimum*** of geometry where there is a difference. 

This is a form of union of the two volumes. 

<p align="center"><img img width="320" src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_min.png?at=refs%2Fheads%2Fmaster"/></p>


