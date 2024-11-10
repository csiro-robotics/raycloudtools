# NetCDF_CXXConfig.cmake
#
# Finds the NetCDF C++ library
#
# This will define the following variables:
#
#    NetCDF_CXX_FOUND        - True if the NetCDF C++ library is found
#    NetCDF_CXX_LIBRARIES    - The NetCDF C++ libraries
#    NetCDF_CXX_INCLUDE_DIRS - The NetCDF C++ include directories
#    NetCDF_CXX_VERSION      - The NetCDF C++ version string
#
# and the following imported target:
#
#     NetCDF::CXX

include(CMakeFindDependencyMacro)

# Try to find NetCDF C library first, as C++ library depends on it
find_dependency(NetCDF REQUIRED)

# Find the NetCDF C++ library
find_library(NetCDF_CXX_LIBRARY
    NAMES netcdf_c++4 netcdf-cxx4
    PATHS $ENV{NetCDF_CXX_ROOT}/lib ${NetCDF_CXX_ROOT}/lib
)

# Find the NetCDF C++ headers
find_path(NetCDF_CXX_INCLUDE_DIR
    NAMES netcdf
    PATHS $ENV{NetCDF_CXX_ROOT}/include ${NetCDF_CXX_ROOT}/include
)

# Extract version information
if(NetCDF_CXX_INCLUDE_DIR)
    file(STRINGS "${NetCDF_CXX_INCLUDE_DIR}/netcdf" netcdf_version_str
         REGEX "^#define[ \t]+NETCDF_VERSION[ \t]+\"[^\"]*\"")
    if(netcdf_version_str)
        string(REGEX REPLACE "^#define[ \t]+NETCDF_VERSION[ \t]+\"([^\"]*)\".*" "\\1"
               NetCDF_CXX_VERSION "${netcdf_version_str}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF_CXX
    REQUIRED_VARS NetCDF_CXX_LIBRARY NetCDF_CXX_INCLUDE_DIR
    VERSION_VAR NetCDF_CXX_VERSION
)

if(NetCDF_CXX_FOUND)
    set(NetCDF_CXX_LIBRARIES ${NetCDF_CXX_LIBRARY})
    set(NetCDF_CXX_INCLUDE_DIRS ${NetCDF_CXX_INCLUDE_DIR})

    if(NOT TARGET NetCDF::CXX)
        add_library(NetCDF::CXX UNKNOWN IMPORTED)
        set_target_properties(NetCDF::CXX PROPERTIES
            IMPORTED_LOCATION "${NetCDF_CXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_CXX_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES NetCDF::NetCDF
        )
    endif()
endif()

mark_as_advanced(NetCDF_CXX_LIBRARY NetCDF_CXX_INCLUDE_DIR)