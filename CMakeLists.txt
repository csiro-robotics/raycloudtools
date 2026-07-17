# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas
# Author: Rod Taylor

cmake_minimum_required(VERSION 3.22) # Ubuntu 22.04 has CMake 3.22.

# Setup project details.
project(raycloudtools VERSION 0.1.0
  DESCRIPTION "A collection of tools for working with point cloud data, including rayclouds."
  LANGUAGES CXX
)

# Add project level options here.
option(RAYCLOUD_BUILD_TESTS "Build unit tests?" OFF)
set(RAYLIB_IMAGETOOL kolourpaint CACHE STRING "command to use to visualise image files")
set(RAYLIB_VISTOOL meshlab CACHE STRING "command to use to visualise point cloud / mesh files")

# WITH_ build options.
option(RAYLIB_DOUBLE_RAYS "Store ray ends as doubles, so distances can be large" OFF)
option(RAYLIB_WITH_LAS "With LASzip for las/laz file support (LAS 1.0-1.4)?" OFF)
option(RAYLIB_WITH_NORMAL_FIELD "Stores rays in the PLY normal field nx,ny,nz. Else rayx,rayy,rayz" ON)
option(RAYLIB_WITH_QHULL "With libqhull support?" OFF)
option(RAYLIB_WITH_TBB "With Intel Threading Building Blocks support multi-threadding?" OFF)
option(RAYLIB_WITH_TIFF "With libgeotiff support?" OFF)

# Setup CMake options.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Setup compiler warnings and options.
add_compile_options(
  "-pedantic"
  "-Wall"
  "-Wextra"
  "-Wconversion"
  "-Werror=pedantic"
  "-Werror=vla"
  "-fdiagnostics-color=always"
)

# Required packages.
find_package(PkgConfig REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

set(RAYTOOLS_INCLUDE
  ${EIGEN3_INCLUDE_DIRS}
  ${libnabo_INCLUDE_DIRS}
)
set(RAYTOOLS_LINK
  libnabo::nabo
  OpenMP::OpenMP_CXX
  Threads::Threads
)

# Optionally configured packages.
if(RAYLIB_WITH_LAS)
  find_package(LASzip REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${LASzip_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${LASzip_LIBRARIES})
endif(RAYLIB_WITH_LAS)

if(RAYLIB_WITH_QHULL)
  find_package(Qhull REQUIRED)
  get_target_property(QHULL_INCLUDE_DIRS Qhull::qhullcpp INTERFACE_INCLUDE_DIRECTORIES)
  set(QHULL_LIBRARIES Qhull::qhullcpp)
  list(APPEND RAYTOOLS_INCLUDE ${QHULL_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${QHULL_LIBRARIES})
endif(RAYLIB_WITH_QHULL)

if(RAYLIB_WITH_TIFF)
  # Find and configure GeoTIFF
  pkg_check_modules(GeoTIFF REQUIRED libgeotiff)
  list(APPEND RAYTOOLS_INCLUDE ${GeoTIFF_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${GeoTIFF_LIBRARIES})

  # Find and configure TIFF
  find_package(TIFF REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${TIFF_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${TIFF_LIBRARIES})

  # Find and configure PROJ
  find_package(PROJ REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${PROJ_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${PROJ_LIBRARIES})
endif(RAYLIB_WITH_TIFF)

if(RAYLIB_WITH_TBB)
  find_package(TBB REQUIRED)
  get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
  set(TBB_LIBRARIES TBB::tbb)
  list(APPEND RAYTOOLS_INCLUDE ${TBB_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${TBB_LIBRARIES})
endif(RAYLIB_WITH_TBB)

# Create libraries and executables.
add_subdirectory(3rd-party)
add_subdirectory(raylib)
add_subdirectory(raycloudtools)

# Export targets.
include(CMakePackageConfigHelpers)
install(
  EXPORT "${CMAKE_PROJECT_NAME}-targets"
  DESTINATION "lib/cmake/${CMAKE_PROJECT_NAME}"
  NAMESPACE "raycloud::"
)
configure_package_config_file(
  "cmake/Config.in.cmake"
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "lib/cmake/${CMAKE_PROJECT_NAME}"
)
install(
  FILES "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-config.cmake"
  DESTINATION "lib/cmake/${CMAKE_PROJECT_NAME}"
)

# Test setup.
if(RAYCLOUD_BUILD_TESTS)
  find_package(GTest REQUIRED)
  enable_testing()
  add_subdirectory(tests)
endif(RAYCLOUD_BUILD_TESTS)
