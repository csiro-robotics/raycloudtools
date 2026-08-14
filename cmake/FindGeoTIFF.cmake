# This module searches for the standalone GeoTIFF library and defines:
#   GeoTIFF_LIBRARIES    - link libraries
#   GeoTIFF_INCLUDE_DIRS - include directories (contains geotiff/geotiff.h)
#   GeoTIFF_FOUND        - true if found

find_path(GeoTIFF_INCLUDE_DIRS
  NAMES geotiff.h
  HINTS ENV GeoTIFF_ROOT
  PATH_SUFFIXES include
)

find_library(GeoTIFF_LIBRARY
  NAMES geotiff
  HINTS ENV GeoTIFF_ROOT
  PATH_SUFFIXES lib
)

set(GeoTIFF_LIBRARIES ${GeoTIFF_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GeoTIFF
  REQUIRED_VARS
    GeoTIFF_INCLUDE_DIRS
    GeoTIFF_LIBRARIES
)

mark_as_advanced(
  GeoTIFF_INCLUDE_DIRS
  GeoTIFF_LIBRARIES
  GeoTIFF_LIBRARY
)
