@PACKAGE_INIT@

# find_dependency.
include(CMakeFindDependencyMacro)

# Let CMake find included Find*.cmake files.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Required dependencies.
find_dependency(Eigen3 REQUIRED)
find_dependency(libnabo REQUIRED)
find_dependency(OpenMP REQUIRED)
find_dependency(Threads REQUIRED)

# Optional dependencies.
set(RAYLIB_WITH_LAS @RAYLIB_WITH_LAS@)
if(RAYLIB_WITH_LAS)
  find_dependency(LASzip REQUIRED)
endif()
set(RAYLIB_WITH_QHULL @RAYLIB_WITH_QHULL@)
if(RAYLIB_WITH_QHULL)
  find_dependency(Qhull REQUIRED)
endif()
set(RAYLIB_WITH_TBB @RAYLIB_WITH_TBB@)
if(RAYLIB_WITH_TBB)
  find_dependency(TBB REQUIRED)
endif()
set(RAYLIB_WITH_TIFF @RAYLIB_WITH_TIFF@)
if(RAYLIB_WITH_TIFF)
  find_dependency(GeoTIFF REQUIRED)
  find_dependency(TIFF REQUIRED)
  find_dependency(PROJ REQUIRED)
endif()

# Targets.
include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
check_required_components("@PROJECT_NAME@")
