#!/usr/bin/env bash
set -o errexit

# Constants.
RAYCLOUDTOOLS_ROOT="$(realpath "$(dirname "${BASH_SOURCE}")/../..")"
PARALLEL_WORKERS="$(("$(nproc)" - 2))"

# Options.
if test -z "${EXAMPLES_ENABLE}"; then
    EXAMPLES_ENABLE=true
fi
if test -z "${RAYCLOUDTOOLS_TEST}"; then
    RAYCLOUDTOOLS_TEST=ON
fi
if test -z "${TREETOOLS_ENABLE}"; then
    if test -d "${RAYCLOUDTOOLS_ROOT}/treetools"; then
        TREETOOLS_ENABLE=true
    else
        TREETOOLS_ENABLE=false
    fi
fi
if test -z "${TREETOOLS_TEST}"; then
    TREETOOLS_TEST=ON
fi

# Check for treetools.
if test "${TREETOOLS_ENABLE}" = "true" && ! test -d "${RAYCLOUDTOOLS_ROOT}/treetools"; then
    echo "TREETOOLS_ENABLE is true but treetools could not be found."
    echo "To use treetools, clone it to ${RAYCLOUDTOOLS_ROOT}/treetools"
    echo "To not use treetools, export TREETOOLS_ENABLE=false"
    echo "https://github.com/csiro-robotics/treetools"
    exit 1
fi

# CMake options.
declare -a CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo"
    "-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON"
    "-DCMAKE_INSTALL_PREFIX:STRING=${RAYCLOUDTOOLS_ROOT}/install"
)
declare -a RAYCLOUDTOOLS_CMAKE_ARGS=(
    "${CMAKE_ARGS[@]}"
    "-DRAYCLOUD_BUILD_TESTS:BOOL=${RAYCLOUDTOOLS_TEST}"
    "-DRAYLIB_DOUBLE_RAYS:BOOL=ON"
    "-DRAYLIB_WITH_LAS:BOOL=ON"
    "-DRAYLIB_WITH_NORMAL_FIELD:BOOL=ON"
    "-DRAYLIB_WITH_QHULL:BOOL=ON"
    "-DRAYLIB_WITH_TBB:BOOL=ON"
    "-DRAYLIB_WITH_TIFF:BOOL=ON"
)
declare -a TREETOOLS_CMAKE_ARGS=(
    "${CMAKE_ARGS[@]}"
    "-DTREE_BUILD_DOXYGEN:BOOL=ON"
    "-DTREE_BUILD_TESTS:BOOL=${TREETOOLS_TEST}"
    "-DWITH_TIFF:BOOL=ON"
)
