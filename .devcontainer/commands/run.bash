#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${RAYCLOUDTOOLS_ROOT}/install/lib"
export PATH="${PATH}:${RAYCLOUDTOOLS_ROOT}/install/bin"

"$@"
