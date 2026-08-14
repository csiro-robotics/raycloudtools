# This module searches for the standalone PROJ library and defines:
#   PROJ_LIBRARIES    - link libraries
#   PROJ_INCLUDE_DIRS - include directories (contains proj.h)
#   PROJ_FOUND        - true if found

find_package(PkgConfig REQUIRED)
pkg_check_modules(PROJ REQUIRED proj)
set(PROJ_INCLUDE_DIRS "${PROJ_INCLUDEDIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PROJ
  REQUIRED_VARS
    PROJ_INCLUDE_DIRS
    PROJ_LIBRARIES
)

mark_as_advanced(
  PROJ_INCLUDE_DIRS
  PROJ_LIBRARIES
)
