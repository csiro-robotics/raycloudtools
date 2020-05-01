# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas
cmake_minimum_required(VERSION 3.10)

# Setup project details.
project(useraylib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(raycloudtools)

add_executable(useraylib useraylib.cpp)
target_link_libraries(useraylib raycloud::raylib)
