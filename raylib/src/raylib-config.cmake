get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/raylib-targets.cmake)
get_filename_component(raylib_INCLUDE_DIRS "${SELF_DIR}/../../../include/raylib" ABSOLUTE)
set_property(TARGET raylib PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${raylib_INCLUDE_DIRS})

set(raylib_FOUND true)
