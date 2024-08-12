#!/bin/bash
set -e

# Print starting message
echo "Starting conditional copying of RiVLib and rdblib folders..."

# Define the base path for Riegl libraries
RIEGL_LIBS_PATH="/workspaces/raycloudtools/riegl_libs"

# Check for the rdblib folder
RDBLIB_PATH=$(find "$RIEGL_LIBS_PATH" -maxdepth 1 -type d -name 'rdblib-*-x86_64-linux' | head -n 1)
if [ -n "$RDBLIB_PATH" ]; then
    echo "Found rdblib folder: $RDBLIB_PATH. Copying to /opt/rdblib..."
    mkdir -p /opt/rdblib
    cp -r "$RDBLIB_PATH"/* /opt/rdblib/
else
    echo "rdblib folder not found."
fi

# Check for the rivlib folder
RIVLIB_PATH=$(find "$RIEGL_LIBS_PATH" -maxdepth 1 -type d -name 'rivlib-*-x86_64-linux-gcc11' | head -n 1)
if [ -n "$RIVLIB_PATH" ]; then
    echo "Found rivlib folder: $RIVLIB_PATH. Copying to /opt/rivlib..."
    mkdir -p /opt/rivlib
    cp -r "$RIVLIB_PATH"/* /opt/rivlib/
else
    echo "rivlib folder not found."
fi

# Print completion message
echo "Conditional copying of RiVLib and rdblib folders completed."
