Ray Cloud Tools

This software is a set of command line tools for processing ray clouds, together with an associated library for using these tools in code.


mkdir build

cd build

cmake ..

make


The command line tools are placed in the bin directory, and the shared library is in the raylib/lib folder.

To access the tools from anywhere, place in your ~/bashrc:
  export PATH=$PATH:<source code path>/raycloudtools/bin
