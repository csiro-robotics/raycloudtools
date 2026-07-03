#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace

# raycloudtools
mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
cmake \
    -S "${RAYCLOUDTOOLS_ROOT}" \
    -B "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
    "${RAYCLOUDTOOLS_CMAKE_ARGS[@]}"
cd "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
ctest \
    --output-junit test_results.gtest.xml \
    --output-on-failure \
    --timeout 300 \
    "$@"
