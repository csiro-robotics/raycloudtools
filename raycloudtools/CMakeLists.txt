# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Rod Taylor

function(add_raycloudtools_executable EXEC_NAME EXEC_SOURCES_LIST)

  add_executable(${EXEC_NAME} ${EXEC_SOURCES_LIST})
  target_link_libraries(${EXEC_NAME} PRIVATE raylib)
  set_target_properties(${EXEC_NAME} PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
  set_target_properties(${EXEC_NAME} PROPERTIES VERSION ${${CMAKE_PROJECT_NAME}_VERSION})
  install(TARGETS ${EXEC_NAME}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION "${INCLUDE_INSTALL_DIR}"
  )

endfunction(add_raycloudtools_executable)


add_subdirectory(rayalign)
add_subdirectory(raycolour)
add_subdirectory(raycombine)
add_subdirectory(raycreate)
add_subdirectory(raydecimate)
add_subdirectory(raydenoise)
add_subdirectory(raydiff)
add_subdirectory(rayexport)
add_subdirectory(rayextract)
add_subdirectory(rayimport)
add_subdirectory(rayinfo)
add_subdirectory(rayrotate)
add_subdirectory(raysmooth)
add_subdirectory(raysplit)
add_subdirectory(raytransients)
add_subdirectory(raytranslate)
add_subdirectory(rayrender)
add_subdirectory(rayrestore)
if(WITH_QHULL)
  add_subdirectory(raywrap)
endif(WITH_QHULL)

