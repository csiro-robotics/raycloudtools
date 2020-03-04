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
