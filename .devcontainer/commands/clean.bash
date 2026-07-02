#!/usr/bin/env bash
set -o errexit
source "$(dirname "${BASH_SOURCE}")/setup.bash"

set -o xtrace
rm --force --recursive \
    "${RAYCLOUDTOOLS_ROOT}/build" \
    "${RAYCLOUDTOOLS_ROOT}/install"
