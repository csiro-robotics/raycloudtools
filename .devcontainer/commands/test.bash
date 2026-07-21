#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace

# raycloudtools.
if test "${RAYCLOUDTOOLS_TEST}" = "ON"; then
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
elif test "${RAYCLOUDTOOLS_TEST}" = "OFF"; then
    echo "Skipping raycloudtools."
else
    echo "RAYCLOUDTOOLS_TEST must be ON or OFF but it is ${RAYCLOUDTOOLS_TEST}."
fi

# treetools.
if test "${TREETOOLS_ENABLE}" = "true"; then
    if test "${TREETOOLS_TEST}" = "ON"; then
        mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/treetools"
        cmake \
            -S "${RAYCLOUDTOOLS_ROOT}/treetools" \
            -B "${RAYCLOUDTOOLS_ROOT}/build/treetools" \
            "${TREETOOLS_CMAKE_ARGS[@]}"
        cd "${RAYCLOUDTOOLS_ROOT}/build/treetools"
        ctest \
            --output-junit test_results.gtest.xml \
            --output-on-failure \
            --timeout 300 \
            "$@"
    elif test "${TREETOOLS_TEST}" = "OFF"; then
        echo "Skipping treetools."
    else
        echo "TREETOOLS_TEST must be ON or OFF but it is ${TREETOOLS_TEST}."
    fi
elif test "${TREETOOLS_ENABLE}" = "false"; then
    echo "Skipping treetools."
else
    echo "TREETOOLS_ENABLE must be true or false but it is ${TREETOOLS_ENABLE}."
fi
