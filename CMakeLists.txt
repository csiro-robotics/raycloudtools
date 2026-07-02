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
# Setup doxygen
option(RAYCLOUD_BUILD_DOXYGEN "Build doxgen documentation?" OFF)
# Setup unit tests
option(RAYCLOUD_BUILD_TESTS "Build unit tests?" OFF)
# Setup LeakTrack
option(RAYCLOUD_LEAK_TRACK "Enable memory leak tracking?" OFF)
set(RAYLIB_VISTOOL meshlab CACHE STRING "command to use to visualise point cloud / mesh files")
set(RAYLIB_IMAGETOOL kolourpaint CACHE STRING "command to use to visualise image files")

# WITH_ build options.
option(RAYLIB_WITH_LAS "With LASzip for las/laz file support (LAS 1.0-1.4)?" OFF)
option(RAYLIB_WITH_QHULL "With libqhull support?" OFF)
option(RAYLIB_WITH_TIFF "With libgeotiff support?" OFF)
option(RAYLIB_WITH_TBB "With Intel Threading Building Blocks support multi-threadding?" OFF)
option(RAYLIB_WITH_NORMAL_FIELD "Stores rays in the PLY normal field nx,ny,nz. Else rayx,rayy,rayz" ON)
option(RAYLIB_DOUBLE_RAYS "Store ray ends as doubles, so distances can be large" OFF)

# Setup CMake options.
set(CMAKE_DEBUG_POSTFIX "d")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(GenerateExportHeader)

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
find_package(Eigen3 REQUIRED)
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads)

set(RAYTOOLS_INCLUDE ${EIGEN3_INCLUDE_DIRS} ${libnabo_INCLUDE_DIRS})
set(RAYTOOLS_LINK ${libnabo_LIBRARIES} OpenMP::OpenMP_CXX Threads::Threads)

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
  find_package(GeoTIFF REQUIRED)

  if(GeoTIFF_FOUND)
    list(APPEND RAYTOOLS_INCLUDE ${GeoTIFF_INCLUDE_DIR})
    list(APPEND RAYTOOLS_LINK ${GeoTIFF_LIBRARIES})
  else()
    message(FATAL_ERROR "GeoTIFF library not found.")
  endif()

  # Find and configure TIFF
  find_package(TIFF REQUIRED)
  if(TIFF_FOUND)
    list(APPEND RAYTOOLS_INCLUDE ${TIFF_INCLUDE_DIRS})
    list(APPEND RAYTOOLS_LINK ${TIFF_LIBRARIES})
  else()
    message(FATAL_ERROR "TIFF library not found.")
  endif()

  # Find and configure PROJ
  find_package(PROJ REQUIRED)
  if(PROJ_FOUND)
    list(APPEND RAYTOOLS_INCLUDE ${PROJ_INCLUDE_DIRS})
    list(APPEND RAYTOOLS_LINK ${PROJ_LIBRARIES})
  else()
    message(FATAL_ERROR "PROJ library not found.")
  endif()
endif(RAYLIB_WITH_TIFF)

if(RAYLIB_WITH_TBB)
  # Helper macro to Find TBB preferring TBBConfig, falling back to FindTBB.
  # TBB_LIBRARIES will be defined in both cases.
  # TBB_INCLUDE_DIRS will be empty for TBBConfig, and populated for FindTBB.
  macro(find_tbb)
    # Try find TBB config file.
    find_package(TBB QUIET CONFIG)
    if(TBB_FOUND)
      message(STATUS "TBB found using TBBConfig.cmake")
      set(TBB_LIBRARIES "${TBB_IMPORTED_TARGETS}")
    else(TBB_FOUND)
      # No TBB config file. Try FindTBB.cmake
      find_package(TBB QUIET)
      if(TBB_FOUND)
        message(STATUS "TBB found using FindTBB.cmake")
      endif(TBB_FOUND)
    endif(TBB_FOUND)
  endmacro(find_tbb)
  # First try newer TBB versions which support TBBConfig.cmake
  find_package(TBB QUIET CONFIG)
  if(TBB_FOUND)
    # Found using TBBConfig. Use import targets
    list(APPEND RAYTOOLS_LINK ${TBB_IMPORTED_TARGETS})
  else(TBB_FOUND)
    # Failed. Fall back to FindTBB.cmake
    find_package(TBB REQUIRED)
    list(APPEND RAYTOOLS_INCLUDE ${TBB_INCLUDE_DIRS})
    list(APPEND RAYTOOLS_LINK ${TBB_LIBRARIES})
  endif(TBB_FOUND)
endif(RAYLIB_WITH_TBB)

# Create libs
add_subdirectory(3rd-party)
add_subdirectory(raylib)
add_subdirectory(raycloudtools)

# Create apps
# add_subdirectory(apps)

# Test setup.
if(RAYCLOUD_BUILD_TESTS)
  find_package(GTest REQUIRED)
  enable_testing()
  add_subdirectory(tests)
endif(RAYCLOUD_BUILD_TESTS)

# Doxygen setup.
if(RAYCLOUD_BUILD_DOXYGEN)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)

    set(DOXYGEN_PROJECT_VERSION ${raycloudtools_VERSION})
    set(DOXYGEN_DOCSET_BUNDLE_ID "au.csiro.raycloudtools")
    set(DOXYGEN_DOCSET_PUBLISHER_ID "au.csiro")
    set(DOXYGEN_DOCSET_PUBLISHER_NAME "CSIRO")

    doxygen_add_docs(
      src
      doc
      ALL
    )
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc/${CMAKE_PROJECT_NAME} COMPONENT Devel)

  else(DOXYGEN_FOUND)
    message(FATAL_ERROR "Told to build Doxygen Documentation, but failed to find Doxygen")
  endif(DOXYGEN_FOUND)

endif(RAYCLOUD_BUILD_DOXYGEN)

# Installation
include(InstallRequiredSystemLibraries)
include(CMakePackageConfigHelpers)

install(EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE ${CMAKE_PROJECT_NAME}-targets.cmake
  NAMESPACE ${PACKAGE_NAMESPACE}
  DESTINATION lib/cmake/${CMAKE_PROJECT_NAME}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-version.cmake"
  VERSION ${raycloudtools_VERSION}
  COMPATIBILITY SameMajorVersion
  #COMPATIBILITY <AnyNewerVersion|SameMajorVersion|SameMinorVersion|ExactVersion>
)

configure_package_config_file(
    "cmake/Config.in.cmake"
    "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-config.cmake"
    INSTALL_DESTINATION "lib/cmake/${CMAKE_PROJECT_NAME}"
)

# Install package files
install(FILES
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-config.cmake"
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-version.cmake"
  DESTINATION "lib/cmake/${CMAKE_PROJECT_NAME}")
