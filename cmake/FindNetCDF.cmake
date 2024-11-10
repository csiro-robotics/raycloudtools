# FindNetCDF.cmake
#
# Finds the NetCDF library
#
# This will define the following variables
#
#    NetCDF_FOUND
#    NetCDF_INCLUDE_DIRS
#    NetCDF_LIBRARIES
#    NetCDF_VERSION
#    NetCDF_HAS_CXX
#
# and the following imported targets
#
#     NetCDF::NetCDF
#     NetCDF::NetCDF_CXX (if C++ library is found)

include(FindPackageHandleStandardArgs)

# Try to find NetCDF with pkg-config (for UNIX-like systems)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_NetCDF QUIET netcdf)
    pkg_check_modules(PC_NetCDF_CXX QUIET netcdf-cxx4)
endif()

# Find the C headers
find_path(NetCDF_INCLUDE_DIR
    NAMES netcdf.h
    PATHS ${PC_NetCDF_INCLUDE_DIRS}
    PATH_SUFFIXES include
)

# Find the C library
find_library(NetCDF_LIBRARY
    NAMES netcdf
    PATHS ${PC_NetCDF_LIBRARY_DIRS}
)

# Find the CXX headers
find_path(NetCDF_CXX_INCLUDE_DIR
    NAMES netcdf
    PATHS ${PC_NetCDF_CXX_INCLUDE_DIRS}
    PATH_SUFFIXES include
)

# Find the CXX library
find_library(NetCDF_CXX_LIBRARY
    NAMES netcdf-cxx4 netcdf_c++4 netcdf-cxx
    PATHS ${PC_NetCDF_CXX_LIBRARY_DIRS}
)

# Extract version information
if(PC_NetCDF_VERSION)
    set(NetCDF_VERSION ${PC_NetCDF_VERSION})
else()
    if(NetCDF_INCLUDE_DIR AND EXISTS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h")
        file(STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" netcdf_version_str
             REGEX "^#define[ \t]+NC_VERSION[ \t]+\"[^\"]*\"")
        string(REGEX REPLACE "^#define[ \t]+NC_VERSION[ \t]+\"([^\"]*)"
               "\\1" NetCDF_VERSION "${netcdf_version_str}")
    endif()
endif()

# Set output variables
set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
set(NetCDF_LIBRARIES ${NetCDF_LIBRARY})

# Check if C++ library is found
if(NetCDF_CXX_LIBRARY AND NetCDF_CXX_INCLUDE_DIR)
    set(NetCDF_HAS_CXX TRUE)
    list(APPEND NetCDF_INCLUDE_DIRS ${NetCDF_CXX_INCLUDE_DIR})
    list(APPEND NetCDF_LIBRARIES ${NetCDF_CXX_LIBRARY})
else()
    set(NetCDF_HAS_CXX FALSE)
endif()

# Handle the QUIETLY and REQUIRED arguments and set NetCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
    VERSION_VAR NetCDF_VERSION
)

# Create imported targets
if(NetCDF_FOUND AND NOT TARGET NetCDF::NetCDF)
    add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF PROPERTIES
        IMPORTED_LOCATION "${NetCDF_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}"
    )
endif()

if(NetCDF_HAS_CXX AND NOT TARGET NetCDF::NetCDF_CXX)
    add_library(NetCDF::NetCDF_CXX UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_CXX PROPERTIES
        IMPORTED_LOCATION "${NetCDF_CXX_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_CXX_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(NetCDF_INCLUDE_DIR NetCDF_CXX_INCLUDE_DIR NetCDF_LIBRARY NetCDF_CXX_LIBRARY)