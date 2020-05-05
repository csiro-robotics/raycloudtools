# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

add_subdirectory(rayalign)
add_subdirectory(raycolour)
add_subdirectory(raycombine)
add_subdirectory(rayconvert)
add_subdirectory(raycreate)
add_subdirectory(raydecimate)
add_subdirectory(raydenoise)
add_subdirectory(rayrotate)
add_subdirectory(raysmooth)
add_subdirectory(raysplit)
add_subdirectory(raytransients)
add_subdirectory(raytranslate)
if(WITH_QHULL)
  add_subdirectory(raywrap)
endif(WITH_QHULL)
