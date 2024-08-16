# FindPROJ.cmake

if (PROJ_FOUND)
    # PROJ is already found, no need to search again
    return()
endif()

# Find the PROJ library
find_package_handle_standard_args(PROJ DEFAULT_MSG PROJ_LIBRARY PROJ_INCLUDE_DIR)

# Define the PROJ library
find_library(PROJ_LIBRARY
    NAMES proj
    PATHS ${CMAKE_LIBRARY_PATH}
    NO_DEFAULT_PATH
)

# Define the PROJ include directory
find_path(PROJ_INCLUDE_DIR
    NAMES proj.h
    PATHS ${CMAKE_INCLUDE_PATH}
    NO_DEFAULT_PATH
)

# Handle the results
include(FindPackageHandleStandardArgs)
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
