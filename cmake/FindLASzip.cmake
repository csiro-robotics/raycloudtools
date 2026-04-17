# This module searches for the standalone LASzip library (3.x+, with LAS 1.4 support)
# and defines:
#   LASzip_LIBRARIES    - link libraries
#   LASzip_INCLUDE_DIRS - include directories (contains laszip/laszip_api.h)
#   LASzip_FOUND        - true if found

find_path(LASzip_INCLUDE_DIRS
  NAMES laszip/laszip_api.h
  HINTS ENV LASzip_ROOT
  PATH_SUFFIXES include
)

find_library(LASzip_LIBRARY
  NAMES laszip
  HINTS ENV LASzip_ROOT
  PATH_SUFFIXES lib
)

set(LASzip_LIBRARIES ${LASzip_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LASzip
  REQUIRED_VARS LASzip_LIBRARIES LASzip_INCLUDE_DIRS
)

mark_as_advanced(LASzip_INCLUDE_DIRS LASzip_LIBRARIES LASzip_LIBRARY)
