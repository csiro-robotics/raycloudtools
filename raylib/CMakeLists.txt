# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# Setup configuration header
configure_file(raylibconfig.in.h "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h")

set(PUBLIC_HEADERS
  rayalignment.h
  rayaxisalign.h
  raycloud.h
  raycloudwriter.h
  rayconcavehull.h
  rayconvexhull.h
  raydecimation.h
  raydiffer.h
  rayellipsoid.h
  rayfinealignment.h
  rayforestgen.h
  rayforeststructure.h
  raygrid.h
  raylaz.h
  raymerger.h
  raymesh.h
  rayply.h
  raypose.h
  rayprogress.h
  rayprogressthread.h
  rayroomgen.h
  raysplitter.h
  raybuildinggen.h
  raycuboid.h
  rayterraingen.h
  raythreads.h
  raytrajectory.h
  raytreegen.h
  raytreestructure.h
  rayunused.h
  rayutils.h
  rayparse.h
  rayrandom.h
  rayrenderer.h
  extraction/raytrunk.h
  extraction/raytrunks.h
  extraction/rayclusters.h
  extraction/rayleaves.h
  extraction/rayterrain.h
  extraction/raytrees.h
  extraction/rayforest.h
  extraction/raysegment.h
  extraction/raytreenode.h
  extraction/raygrid2d.h
)

set(PRIVATE_HEADERS
  imageread.h
  imagewrite.h
)

set(SOURCES
  ${PUBLIC_HEADERS}
  ${PRIVATE_HEADERS}
  rayalignment.cpp
  rayaxisalign.cpp
  raycloud.cpp
  raycloudwriter.cpp
  rayconcavehull.cpp
  rayconvexhull.cpp
  raydecimation.cpp
  raydiffer.cpp
  rayellipsoid.cpp
  rayfinealignment.cpp
  rayforestgen.cpp
  rayforeststructure.cpp
  raylaz.cpp
  raymerger.cpp
  raymesh.cpp
  rayply.cpp
  rayprogressthread.cpp
  rayroomgen.cpp
  raysplitter.cpp
  raybuildinggen.cpp
  raycuboid.cpp
  rayterraingen.cpp
  raythreads.cpp
  raytrajectory.cpp
  raytreegen.cpp
  raytreestructure.cpp
  rayparse.cpp
  rayrandom.cpp
  rayrenderer.cpp
  extraction/raygrid2d.cpp
  extraction/raytrunk.cpp
  extraction/raytrunks.cpp
  extraction/rayclusters.cpp
  extraction/rayleaves.cpp
  extraction/rayterrain.cpp
  extraction/raytrees.cpp
  extraction/rayforest.cpp
  extraction/rayforest_draw.cpp
  extraction/rayforest_watershed.cpp
  extraction/raysegment.cpp
)

get_target_property(SIMPLE_FFT_INCLUDE_DIRS simple_fft INTERFACE_INCLUDE_DIRECTORIES)
add_compile_options(-fPIC)

if(RAYLIB_WITH_QHULL)
  set(QHULL_LIBS
    Qhull::qhullcpp
    Qhull::qhullstatic_r
  )
else(RAYLIB_WITH_QHULL)
  set(QHULL_LIBS)
endif(RAYLIB_WITH_QHULL)

add_library(raylib SHARED ${SOURCES})

generate_export_header(raylib)

target_compile_options(raylib PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(raylib PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  ${RAYTOOLS_INCLUDE}
  ${SIMPLE_FFT_INCLUDE_DIRS}
)
target_link_libraries(raylib PUBLIC ${RAYTOOLS_LINK} PRIVATE ${QHULL_LIBS})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h DESTINATION include/raylib)

install(FILES ${PUBLIC_HEADERS} DESTINATION include/raylib)

install(TARGETS raylib
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include/raylib
)
install(
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  DESTINATION lib/${CMAKE_PROJECT_NAME}
  NAMESPACE raycloud::
)

# ras_add_library(raylib
#   TYPE SHARED
#   INCLUDE_PREFIX "raylib"
#   PROJECT_FOLDER "raylib"
#   INCLUDE
#     PUBLIC_SYSTEM
#       ${RAYTOOLS_INCLUDE}
#     PRIVATE_SYSTEM
#       # Add the simple_fft include directories as PRIVATE system headers to prevent them generating compiler warnings.
#       # The simple_fft dependency is also added as an INTERFACE library, which will be propagated in the dependency
#       # chain.
#       "${SIMPLE_FFT_INCLUDE_DIRS}"
#   LIBS
#     PUBLIC
#       ${RAYTOOLS_LINK}
#     PRIVATE
#       ${QHULL_LIBS}
#   PUBLIC_HEADERS ${PUBLIC_HEADERS}
#   GENERATED PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h"
#   SOURCES ${SOURCES}
# )
