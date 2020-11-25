# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

# Setup configuration header
configure_file(raylibconfig.in.h "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h")

set(PUBLIC_HEADERS
  rayalignment.h
  raycloud.h
  raycloudwriter.h
  rayconcavehull.h
  rayconvexhull.h
  raydebugdraw.h
  rayellipsoid.h
  rayfinealignment.h
  rayforestgen.h
  raygrid.h
  raylaz.h
  raymerger.h
  raymesh.h
  rayply.h
  raypose.h
  rayprogress.h
  rayprogressthread.h
  rayroomgen.h
  raybuildinggen.h
  raycuboid.h
  rayterraingen.h
  raythreads.h
  raytrajectory.h
  raytreegen.h
  rayunused.h
  rayutils.h
  rayparse.h
  rayrandom.h
)

set(PRIVATE_HEADERS
  imagewrite.h
)

set(SOURCES
  ${PUBLIC_HEADERS}
  ${PRIVATE_HEADERS}
  rayalignment.cpp
  raycloud.cpp
  raycloudwriter.cpp
  rayconcavehull.cpp
  rayconvexhull.cpp
  rayellipsoid.cpp
  rayfinealignment.cpp
  rayforestgen.cpp
  raylaz.cpp
  raymerger.cpp
  raymesh.cpp
  rayply.cpp
  rayprogressthread.cpp
  rayroomgen.cpp
  raybuildinggen.cpp
  raycuboid.cpp
  rayterraingen.cpp
  raythreads.cpp
  raytrajectory.cpp
  raytreegen.cpp
  rayparse.cpp
  rayrandom.cpp
)

# Select the source file to use with raydebudraw.
if(WITH_3ES)
  # Using 3rd Eye Scene
  list(APPEND SOURCES raydebugdraw_3es.cpp)
elseif(WITH_ROS)
  # Using ROS/rivz
  list(APPEND SOURCES raydebugdraw_ros.cpp)
else(WITH_3ES)
  # Disabled.
  list(APPEND SOURCES raydebugdraw_none.cpp)
endif(WITH_3ES)

get_target_property(SIMPLE_FFT_INCLUDE_DIRS simple_fft INTERFACE_INCLUDE_DIRECTORIES)

if(WITH_QHULL)
set(QHULL_LIBS
    Qhull::qhullcpp
    Qhull::qhullstatic_r)
else(WITH_QHULL)
set(QHULL_LIBS)
endif(WITH_QHULL)

ras_add_library(raylib
  TYPE SHARED
  INCLUDE_PREFIX "raylib"
  PROJECT_FOLDER "raylib"
  INCLUDE
    PUBLIC_SYSTEM
      ${RAYTOOLS_INCLUDE}
    PRIVATE_SYSTEM
      # Add the simple_fft include directories as PRIVATE system headers to prevent them generating compiler warnings.
      # The simple_fft dependency is also added as an INTERFACE library, which will be propagated in the dependency
      # chain.
      "${SIMPLE_FFT_INCLUDE_DIRS}"
  LIBS
    PUBLIC
      ${RAYTOOLS_LINK}
    PRIVATE
      ${QHULL_LIBS}
  PUBLIC_HEADERS ${PUBLIC_HEADERS}
  GENERATED PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h"
  SOURCES ${SOURCES}
)

target_compile_options(raylib PUBLIC ${OpenMP_CXX_FLAGS})
