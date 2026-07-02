#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace

mkdir --parents "${RAYCLOUDTOOLS_ROOT}/install"

# raycloudtools
mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
cmake \
    -S "${RAYCLOUDTOOLS_ROOT}" \
    -B "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
    "${RAYCLOUDTOOLS_CMAKE_ARGS[@]}"
make \
    -C "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
    -j "${PARALLEL_WORKERS}" \
    install

export CMAKE_MODULE_PATH="${RAYCLOUDTOOLS_ROOT}/install/lib/cmake/raycloudtools"

# treetools
mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/treetools"
cmake \
    -S "${RAYCLOUDTOOLS_ROOT}/treetools" \
    -B "${RAYCLOUDTOOLS_ROOT}/build/treetools" \
    "${TREETOOLS_CMAKE_ARGS[@]}"
make \
    -C "${RAYCLOUDTOOLS_ROOT}/build/treetools" \
    -j "${PARALLEL_WORKERS}" \
    install
