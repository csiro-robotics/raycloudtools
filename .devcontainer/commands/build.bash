#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace

mkdir --parents "${RAYCLOUDTOOLS_ROOT}/install"

# raycloudtools.
mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
cmake \
    -S "${RAYCLOUDTOOLS_ROOT}" \
    -B "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
    "${RAYCLOUDTOOLS_CMAKE_ARGS[@]}"
make \
    -C "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
    -j "${PARALLEL_WORKERS}" \
    install

export CMAKE_MODULE_PATH="${RAYCLOUDTOOLS_ROOT}/install/lib/raycloudtools"

# Examples.
if test "${EXAMPLES_ENABLE}" = "true"; then
    for EXAMPLE_ROOT in "${RAYCLOUDTOOLS_ROOT}/examples/"*; do
        EXAMPLE_NAME="$(basename "${EXAMPLE_ROOT}")"
        echo "${EXAMPLE_NAME}"
        mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/examples/${EXAMPLE_NAME}"
        cmake \
            -S "${RAYCLOUDTOOLS_ROOT}/examples/${EXAMPLE_NAME}" \
            -B "${RAYCLOUDTOOLS_ROOT}/build/examples/${EXAMPLE_NAME}" \
            "${CMAKE_ARGS[@]}"
        make \
            -C "${RAYCLOUDTOOLS_ROOT}/build/examples/${EXAMPLE_NAME}" \
            -j "${PARALLEL_WORKERS}" \
            install
    done
elif test "${EXAMPLES_ENABLE}" = "false"; then
    echo "Skipping examples."
else
    echo "EXAMPLES_ENABLE must be true or false but it is ${EXAMPLES_ENABLE}."
fi

# treetools.
if test "${TREETOOLS_ENABLE}" = "true"; then
    mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/treetools"
    cmake \
        -S "${RAYCLOUDTOOLS_ROOT}/treetools" \
        -B "${RAYCLOUDTOOLS_ROOT}/build/treetools" \
        "${TREETOOLS_CMAKE_ARGS[@]}"
    make \
        -C "${RAYCLOUDTOOLS_ROOT}/build/treetools" \
        -j "${PARALLEL_WORKERS}" \
        install
elif test "${TREETOOLS_ENABLE}" = "false"; then
    echo "Skipping treetools."
else
    echo "TREETOOLS_ENABLE must be true or false but it is ${TREETOOLS_ENABLE}."
fi
