# FindPROJ.cmake

if (PROJ_FOUND)
    # PROJ is already found, no need to search again
    return()
endif()

include(FindPackageHandleStandardArgs)

# Define the PROJ library
find_library(PROJ_LIBRARY
    NAMES proj
        HINTS
            ${CMAKE_INSTALL_PREFIX}/lib
            /usr/local/lib
            /usr/lib
        PATH_SUFFIXES
            ${CMAKE_LIBRARY_ARCHITECTURE}
            aarch64-linux-gnu
)

# Define the PROJ include directory
find_path(PROJ_INCLUDE_DIR
    NAMES proj.h
        HINTS
            ${CMAKE_INSTALL_PREFIX}/include
            /usr/local/include
            /usr/include
        PATH_SUFFIXES
            proj
)

# Handle the results
find_package_handle_standard_args(PROJ
    REQUIRED_VARS PROJ_LIBRARY PROJ_INCLUDE_DIR
    VERSION_VAR PROJ_VERSION
)

# Provide variables for the library and include directory
set(PROJ_INCLUDE_DIRS ${PROJ_INCLUDE_DIR})
set(PROJ_LIBRARIES ${PROJ_LIBRARY})

# Provide version information if available
if (PROJ_VERSION)
    set(PROJ_VERSION_STRING "${PROJ_VERSION}")
else()
    set(PROJ_VERSION_STRING "unknown")
endif()
