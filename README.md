Ray Cloud Tools

This software is a set of command line tools for processing ray clouds, together with an associated library for using these tools in code.


mkdir build

cd build

cmake ..

make


The command line tools are placed in the bin directory, and the shared library is in the raylib/lib folder.

To access the tools from anywhere, place in your ~/bashrc:
  export PATH=$PATH:<source code path>/raycloudtools/bin


To visualise rayclouds use meshlab. To view the rays choose Render menu then Show Vertex Normals. In the Tools menu select options and change the Decoration::NormalLength to 0.0025 to render them approximately the correct length. 

Usage:

<img
img width="640"
src="https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room1.png?at=refs%2Fheads%2Fmaster"
alt="Subject Pronouns"
style="margin-right: 10px;"
/>


![room1](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room1.png?at=refs%2Fheads%2Fmaster)

![room2](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room2.png?at=refs%2Fheads%2Fmaster)

![roominside](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room3.png?at=refs%2Fheads%2Fmaster)

![roomcombinedall](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_all.png?at=refs%2Fheads%2Fmaster)

![roomcombinedmin](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_combined_min.png?at=refs%2Fheads%2Fmaster)

![roomdecimated](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_decimated.png?at=refs%2Fheads%2Fmaster)

![roomdenoise1](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise1.png?at=refs%2Fheads%2Fmaster)

![roomdenoise2](https://bitbucket.csiro.au/projects/ASR/repos/raycloudtools/raw/pics/room_denoise2.png?at=refs%2Fheads%2Fmaster)
