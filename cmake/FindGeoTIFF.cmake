# FindGeoTIFF.cmake

# Try to find the GeoTIFF library and headers
find_path(GeoTIFF_INCLUDE_DIR NAMES geotiff.h
    HINTS
    ${CMAKE_INSTALL_PREFIX}/include
    /usr/local/include/geotiff
    /usr/include
)

find_library(GeoTIFF_LIBRARY NAMES geotiff
    HINTS
    ${CMAKE_INSTALL_PREFIX}/lib
    /usr/local/lib
    /usr/lib
    /usr/lib/x86_64-linux-gnu
)

# Handle the result
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GeoTIFF REQUIRED_VARS GeoTIFF_INCLUDE_DIR GeoTIFF_LIBRARY)

if(GeoTIFF_FOUND)
    set(GeoTIFF_LIBRARIES ${GeoTIFF_LIBRARY})
else()
    set(GeoTIFF_LIBRARIES "")
endif()

mark_as_advanced(GeoTIFF_INCLUDE_DIR GeoTIFF_LIBRARY)
