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
    "RAYCLOUD_BUILD_DOXYGEN"
    "RAYCLOUD_BUILD_TESTS"
    "RAYCLOUD_LEAK_TRACK"
    "RAYLIB_WITH_LAS"
    "RAYLIB_WITH_QHULL"
    "RAYLIB_WITH_TIFF"
    "RAYLIB_WITH_TBB"
    "RAYLIB_WITH_NORMAL_FIELD"
    "RAYLIB_DOUBLE_RAYS"
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
    rm --force --recursive \
        "${RAYCLOUDTOOLS_ROOT}/build" \
        "${RAYCLOUDTOOLS_ROOT}/install"
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
