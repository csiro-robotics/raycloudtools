#!/bin/bash
set -e

# Build RayCloudTools
echo "Building RayCloudTools..."
cd /workspaces/raycloudtools
mkdir -p build
cd build 

# Construct the cmake command
cmake .. \
    -DGeoTIFF_INCLUDE_DIR=/usr/include/geotiff \
    -DGeoTIFF_LIBRARY=/usr/lib/x86_64-linux-gnu/libgeotiff.so \
    -DPROJ_INCLUDE_DIR=/usr/include/proj \
    -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
    -DWITH_QHULL=ON \
    -DWITH_LAS=ON \
    -DDOUBLE_RAYS=ON \
    -DWITH_TIFF=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DRAYCLOUD_BUILD_TESTS=ON

make -j"$(nproc)"
make install
cd /workspaces/raycloudtools
ldconfig /usr/local/lib

# Clone and build TreeTools
echo "Cloning and building TreeTools..."
rm -rf treetools && git clone https://github.com/csiro-robotics/treetools.git
cd treetools
mkdir -p build
cd build
cmake .. \
    -DPROJ_INCLUDE_DIR=/usr/include/proj \
    -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
    -DWITH_TIFF=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DTREE_BUILD_TESTS=ON
make -j"$(nproc)"
make install
ldconfig /usr/local/lib

echo "Build process completed successfully."