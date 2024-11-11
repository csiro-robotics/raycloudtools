#!/bin/bash
set -e

# Print starting message
echo "Starting conditional copying of RiVLib and rdblib folders..."

# Define the base path for Riegl libraries
RIEGL_LIBS_PATH="/workspaces/raycloudtools/riegl_libs"

# Initialize flags
WITH_RIEGL=OFF

# Check for the rdblib folder
RDBLIB_PATH=$(find "$RIEGL_LIBS_PATH" -maxdepth 1 -type d -name 'rdblib-*-x86_64-linux' | head -n 1)
if [ -n "$RDBLIB_PATH" ]; then
    echo "Found rdblib folder: $RDBLIB_PATH. Copying to /opt/rdblib..."
    mkdir -p /opt/rdblib
    cp -r "$RDBLIB_PATH"/* /opt/rdblib/
    WITH_RIEGL=ON
    # Set environment variables for Riegl libraries
    export RiVLib_DIR=/opt/rivlib
    export RiVLib_INCLUDE_DIRS=/opt/rivlib/include 
else
    echo "rdblib folder not found."
fi

# Check for the rivlib folder
RIVLIB_PATH=$(find "$RIEGL_LIBS_PATH" -maxdepth 1 -type d -name 'rivlib-*-x86_64-linux-gcc11' | head -n 1)
if [ -n "$RIVLIB_PATH" ]; then
    echo "Found rivlib folder: $RIVLIB_PATH. Copying to /opt/rivlib..."
    mkdir -p /opt/rivlib
    cp -r "$RIVLIB_PATH"/* /opt/rivlib/
    WITH_RIEGL=ON
    # Set environment variables for Riegl libraries
    export rdb_DIR=/opt/rdblib/interface/cpp 
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/rdblib/library
else
    echo "rivlib folder not found."
fi

# Print completion message
echo "Conditional copying of RiVLib and rdblib folders completed."

# Build RayCloudTools
echo "Building RayCloudTools..."
cd /workspaces/raycloudtools
mkdir -p build
cd build 

# Construct the cmake command with conditional Riegl support
cmake .. \
    -DGeoTIFF_INCLUDE_DIR=/usr/include/geotiff \
    -DGeoTIFF_LIBRARY=/usr/lib/x86_64-linux-gnu/libgeotiff.so \
    -DPROJ_INCLUDE_DIR=/usr/include/proj \
    -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
    -DWITH_LAS=ON \
    -DDOUBLE_RAYS=ON \
    -DWITH_TIFF=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DWITH_RIEGL=$WITH_RIEGL \
    -DWITH_NETCDF=ON

make -j"$(nproc)"
make install
cd /workspaces/raycloudtools
ldconfig /usr/local/lib

# Clone and build TreeTools
echo "Cloning and building TreeTools..."
git clone https://github.com/Leaf2Landscape/treetools.git
cd treetools
mkdir -p build
cd build
cmake .. \
    -DPROJ_INCLUDE_DIR=/usr/include/proj \
    -DPROJ_LIBRARY=/usr/lib/x86_64-linux-gnu/libproj.so \
    -DWITH_TIFF=ON \
    -DWITH_NETCDF=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j"$(nproc)"
make install
ldconfig /usr/local/lib

echo "Build process completed successfully."
