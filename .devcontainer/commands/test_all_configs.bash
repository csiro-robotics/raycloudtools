#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace

# Print current target on failure.
PROJECT=NONE
ENABLE_OPTION=NONE
function on_exit() {
    echo "PROJECT=${PROJECT}"
    echo "ENABLE_OPTION=${ENABLE_OPTION}"
}
trap on_exit EXIT

# raycloudtools
PROJECT=raycloudtools
declare -a OPTIONS=(
    "RAYLIB_DOUBLE_RAYS"
    "RAYLIB_WITH_LAS"
    "RAYLIB_WITH_NORMAL_FIELD"
    "RAYLIB_WITH_QHULL"
    "RAYLIB_WITH_TBB"
    "RAYLIB_WITH_TIFF"
)
for ENABLE_OPTION in "${OPTIONS[@]}"; do
    declare -a OPTION_CMAKE_ARGS=()
    for OPTION in "${OPTIONS[@]}"; do
        if test "${OPTION}" = "${ENABLE_OPTION}"; then
            OPTION_CMAKE_ARGS+=("-D${OPTION}:BOOL=ON")
        else
            OPTION_CMAKE_ARGS+=("-D${OPTION}:BOOL=OFF")
        fi
    done
    cd "${RAYCLOUDTOOLS_ROOT}"

    # Clean.
    rm --force --recursive \
        "${RAYCLOUDTOOLS_ROOT}/build" \
        "${RAYCLOUDTOOLS_ROOT}/install"

    # raycloudtools.
    mkdir --parents "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
    cmake \
        -S "${RAYCLOUDTOOLS_ROOT}" \
        -B "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
        "${RAYCLOUDTOOLS_CMAKE_ARGS[@]}" \
        "${OPTION_CMAKE_ARGS[@]}"
    make \
        -C "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools" \
        -j "${PARALLEL_WORKERS}" \
        install

    # Make sure options are set.
    for OPTION in "${OPTIONS[@]}"; do
        if test "${OPTION}" = "${ENABLE_OPTION}"; then
            grep "${OPTION}:BOOL=ON" "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools/CMakeCache.txt"
            grep "#define ${OPTION} 1" "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools/raylib/raylibconfig.h"
        else
            grep "${OPTION}:BOOL=OFF" "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools/CMakeCache.txt"
            grep "#define ${OPTION} 0" "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools/raylib/raylibconfig.h"
        fi
    done

    # Examples.
    for EXAMPLE_ROOT in "${RAYCLOUDTOOLS_ROOT}/examples/"*; do
        EXAMPLE_NAME="$(basename "${EXAMPLE_ROOT}")"
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

    # Test.
    mkdir --parents "${RAYCLOUDTOOLS_ROOT}/test_results"
    cd "${RAYCLOUDTOOLS_ROOT}/build/raycloudtools"
    ctest \
        --output-junit "${RAYCLOUDTOOLS_ROOT}/test_results/test_results.${ENABLE_OPTION}.gtest.xml" \
        --output-on-failure \
        --timeout 300
done

# Done.
PROJECT=DONE
ENABLE_OPTION=DONE
