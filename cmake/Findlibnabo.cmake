# - Try to find the libnabo library
# Once done, this will define
#
#  libnabo_FOUND - system has libnabo
#  libnabo_INCLUDE_DIRS - the libnabo include directories
#  libnabo_LIBRARIES - link these to use libnabo
#  libnabo_VERSION - the version of libnabo found

find_path(libnabo_INCLUDE_DIRS
  NAMES nabo/nabo.h
  PATHS
    ${CMAKE_INSTALL_PREFIX}/include
    /usr/local/include
    /usr/include
    /usr/local
    # Add other paths if needed
)

find_library(libnabo_LIBRARIES
  NAMES nabo
  PATHS
    ${CMAKE_INSTALL_PREFIX}/lib
    /usr/local/lib
    /usr/lib
    # Add other paths if needed
)

# Optional: find_package for dependencies, like Eigen
find_package(Eigen3 QUIET)

if(libnabo_INCLUDE_DIRS AND libnabo_LIBRARIES)
  set(libnabo_FOUND TRUE)
else()
  set(libnabo_FOUND FALSE)
endif()

# Provide an imported target if found
if(libnabo_FOUND)
  add_library(libnabo INTERFACE IMPORTED)
  set_target_properties(libnabo PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${libnabo_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${libnabo_LIBRARIES}"
  )
endif()

mark_as_advanced(libnabo_INCLUDE_DIRS libnabo_LIBRARIES)
