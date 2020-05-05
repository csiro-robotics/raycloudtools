# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# We setup simple_fft as an interface library intended only for private, internal usage. It is not installed
# thus not publically visible and must only be included in private header and source files without raycloud tools.
add_library(simple_fft INTERFACE)
# Setup include directories.
target_include_directories(simple_fft INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
