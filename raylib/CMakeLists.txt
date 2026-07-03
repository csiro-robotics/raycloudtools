# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# Setup configuration header
configure_file(raylibconfig.in.h "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h")

set(PUBLIC_HEADERS
  extraction/rayclusters.h
  extraction/rayforest.h
  extraction/raygrid2d.h
  extraction/rayleaves.h
  extraction/raysegment.h
  extraction/rayterrain.h
  extraction/raytreenode.h
  extraction/raytrees.h
  extraction/raytrunk.h
  extraction/raytrunks.h
  rayalignment.h
  rayaxisalign.h
  raybuildinggen.h
  raycloud.h
  raycloudwriter.h
  rayconcavehull.h
  rayconvexhull.h
  raycuboid.h
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
  rayparse.h
  rayply.h
  raypose.h
  rayprogress.h
  rayprogressthread.h
  rayrandom.h
  rayrenderer.h
  rayroomgen.h
  raysplitter.h
  rayterraingen.h
  raythreads.h
  raytrajectory.h
  raytreegen.h
  raytreestructure.h
  rayunused.h
  rayutils.h
)

set(PRIVATE_HEADERS
  imageread.h
  imagewrite.h
)

set(SOURCES
  ${PRIVATE_HEADERS}
  ${PUBLIC_HEADERS}
  extraction/rayclusters.cpp
  extraction/rayforest_draw.cpp
  extraction/rayforest_watershed.cpp
  extraction/rayforest.cpp
  extraction/raygrid2d.cpp
  extraction/rayleaves.cpp
  extraction/raysegment.cpp
  extraction/rayterrain.cpp
  extraction/raytrees.cpp
  extraction/raytrunk.cpp
  extraction/raytrunks.cpp
  rayalignment.cpp
  rayaxisalign.cpp
  raybuildinggen.cpp
  raycloud.cpp
  raycloudwriter.cpp
  rayconcavehull.cpp
  rayconvexhull.cpp
  raycuboid.cpp
  raydecimation.cpp
  raydiffer.cpp
  rayellipsoid.cpp
  rayfinealignment.cpp
  rayforestgen.cpp
  rayforeststructure.cpp
  raylaz.cpp
  raymerger.cpp
  raymesh.cpp
  rayparse.cpp
  rayply.cpp
  rayprogressthread.cpp
  rayrandom.cpp
  rayrenderer.cpp
  rayroomgen.cpp
  raysplitter.cpp
  rayterraingen.cpp
  raythreads.cpp
  raytrajectory.cpp
  raytreegen.cpp
  raytreestructure.cpp
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

include(GenerateExportHeader)
generate_export_header(raylib)

target_compile_options(raylib PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(raylib PUBLIC
  ${RAYTOOLS_INCLUDE}
  ${SIMPLE_FFT_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(raylib PUBLIC ${RAYTOOLS_LINK} PRIVATE ${QHULL_LIBS})

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/raylib_export.h"
    "${CMAKE_CURRENT_BINARY_DIR}/raylibconfig.h"
  DESTINATION include/raylib)

install(FILES ${PUBLIC_HEADERS} DESTINATION include/raylib)

install(TARGETS raylib
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include/raylib
)
