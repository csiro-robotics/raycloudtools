# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas
cmake_minimum_required(VERSION 3.10)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "Build type empty, so defaulting to Release.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Setup project details.
project(raycloudtools)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(RasProjectSetup)
ras_project(
  VERSION 0.0.1
  NAMESPACE raycloud
)

# Add project level options here.
# Setup doxygen
option(RAYCLOUD_BUILD_DOXYGEN "Build doxgen documentation?" OFF)
# Setup unit tests
option(RAYCLOUD_BUILD_TESTS "Build unit tests?" OFF)
# Setup LeakTrack
option(RAYCLOUD_LEAK_TRACK "Enable memory leak tracking?" OFF)

# WITH_ build options.
option(WITH_LAS "With liblas for las file support?" OFF)
option(WITH_QHULL "With libqhull support?" OFF)
option(WITH_TIFF "With libgeotiff support?" OFF)
option(WITH_TBB "With Intel Threading Building Blocks support multi-threadding?" OFF)
option(WITH_NORMAL_FIELD "Stores rays in the PLY normal field nx,ny,nz. Else rayx,rayy,rayz" ON)

# Convert WITH_ options to 1/0 so we can use them in configuration headers.
ras_bool_to_int(WITH_LAS)
ras_bool_to_int(WITH_QHULL)
ras_bool_to_int(WITH_TIFF)
ras_bool_to_int(WITH_TBB)
ras_bool_to_int(WITH_NORMAL_FIELD)

# other build-time options
option(DOUBLE_RAYS "Store ray ends as doubles, so distances can be large" OFF)
ras_bool_to_int(DOUBLE_RAYS)

# Required packages.
find_package(Eigen3 REQUIRED)
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads)

set(RAYTOOLS_INCLUDE ${EIGEN3_INCLUDE_DIRS} ${libnabo_INCLUDE_DIRS})
set(RAYTOOLS_LINK ${libnabo_LIBRARIES} Threads::Threads)

# Optionally configured packages.
if(WITH_LAS)
  find_package(libLAS REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${libLAS_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${libLAS_LIBRARIES})
endif(WITH_LAS)
if(WITH_QHULL)
  # Find the qhull that provides a config with targets first
  find_package(Qhull CONFIG QUIET)
  # Otherwise we try and use our own findQhull.cmake script
  if(NOT Qhull_FOUND)
    find_package(Qhull REQUIRED)
  endif(NOT Qhull_FOUND)
  
  list(APPEND RAYTOOLS_INCLUDE ${QHULL_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${QHULL_LIBRARIES})
endif(WITH_QHULL)
if(WITH_TIFF)
  find_package(GeoTIFF REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${GEOTIFF_INCLUDE_DIR})
  list(APPEND RAYTOOLS_LINK ${GEOTIFF_LIBRARIES})
  find_package(TIFF REQUIRED)
  list(APPEND RAYTOOLS_INCLUDE ${TIFF_INCLUDE_DIRS})
  list(APPEND RAYTOOLS_LINK ${TIFF_LIBRARIES})
  find_package(PROJ REQUIRED CONFIG)
  list(APPEND RAYTOOLS_INCLUDE ${PROJ_INCLUDE_DIR})
  list(APPEND RAYTOOLS_LINK ${PROJ_LIBRARIES})
endif(WITH_TIFF)

if(WITH_TBB)
# Helper macro to Find TBB prefering TBBConfig, falling back to FindTBB. TBB_LIBRARIES will be defined in both cases.
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
endif(WITH_TBB)

# Create libs
add_subdirectory(3rd-party)
add_subdirectory(raylib)
add_subdirectory(raycloudtools)

# Create apps
# add_subdirectory(apps)

# Test setup.
if(RAYCLOUD_BUILD_TESTS)
  find_package(GTest REQUIRED)

  # We can enable testing here and/or in the subdirectory, but doing it here allows us to run CTest from the build root.
  # To run the tests, we execute:
  #   CTest -C [Debug|Release|RelWithDebInfo|MinSizeRel] --output-on-failure
  # CTest normally shows only a very terse test ouput, but we make sure failed tests show all output by adding
  #   --output-on-failure
  # The full test output is always available in:
  #   <build>/Testing/Temporary/LastTest.log
  enable_testing()
  add_subdirectory(tests)
endif(RAYCLOUD_BUILD_TESTS)

# Doxygen setup.
if(RAYCLOUD_BUILD_DOXYGEN)
  # Include Doxygen helper functions. This also finds the Doxygen package.
  include(RasDoxygen)

  if(DOXYGEN_FOUND)
    # Create a target to build the documentation.
    # Here we also setup various documentation variables passed through to the doxyfile configuration.
    # Each named argument below describes the Doxygen variable it sets.
    ras_doxygen_create(
      # DOXYFILE cmake/doxyfile.in  # Doxyfile to configure.
      PROJECT ${CMAKE_PROJECT_NAME} # PROJECT_NAME
      VERSION ${raycloudtools_VERSION}   # PROJECT_NUMBER
      OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html # HTML_OUPTUT
      # CSS <style>.css             # HTML_STYLESHEET
      PUBLISHER "CSIRO"             # DOCSET_PUBLISHER_NAME
      PUBLISHER_ID au.csiro         # DOCSET_PUBLISHER_ID
      PROJECT_ID au.csiro.ras  # DOCSET_BUNDLE_ID, QHP_NAMESPACE, ECLIPSE_DOC_ID
      PATHS                         # INPUT (RECURSIVE is on)
        src
        doc
      EXCLUDE_PATHS                 # EXCLUDE
      # Where to find source code examples.
      # EXAMPLE_PATHS <paths>        # EXAMPLE_PATH
      # Where to find images.
      # IMAGE_PATHS <paths>          # IMAGE_PATH
    )

    # Setup installation of the generated documentation: source, destination.
    ras_doxygen_install(${CMAKE_CURRENT_BINARY_DIR}/html ras)
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
  DESTINATION ${PACKAGE_EXPORT_LOCATION}
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
    INSTALL_DESTINATION "${PACKAGE_EXPORT_LOCATION}"
)

# Install package files
install(FILES
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-config.cmake"
  "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-version.cmake"
  DESTINATION "${PACKAGE_EXPORT_LOCATION}")
