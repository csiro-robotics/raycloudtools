# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# Setup configuration header
configure_file(raylibconfig.in.h "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h")

set(PUBLIC_HEADERS
  rayalignment.h
  raycloud.h
  rayconcavehull.h
  rayconvexhull.h
  raydebugdraw.h
  rayforestgen.h
  raygrid.h
  raylaz.h
  raymesh.h
  rayply.h
  raypose.h
  rayroomgen.h
  rayterraingen.h
  raytrajectory.h
  raytreegen.h
  rayunused.h
  rayutils.h
)

set(PRIVATE_HEADERS
  imagewrite.h
)

set(SOURCES
  ${PUBLIC_HEADERS}
  ${PRIVATE_HEADERS}
  rayalignment.cpp
  raycloud.cpp
  rayconcavehull.cpp
  rayconvexhull.cpp
  raydebugdraw.cpp
  rayforestgen.cpp
  raylaz.cpp
  raymesh.cpp
  rayply.cpp
  rayroomgen.cpp
  rayterraingen.cpp
  raytrajectory.cpp
  raytreegen.cpp
)

get_target_property(SIMPLE_FFT_INCLUDE_DIRS simple_fft INTERFACE_INCLUDE_DIRECTORIES)

ras_add_library(raylib
  TYPE STATIC
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
  LIBS PUBLIC ${RAYTOOLS_LINK}
  PUBLIC_HEADERS ${PUBLIC_HEADERS}
  GENERATED PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h"
  SOURCES ${SOURCES}
)

target_compile_options(raylib PUBLIC ${OpenMP_CXX_FLAGS})
