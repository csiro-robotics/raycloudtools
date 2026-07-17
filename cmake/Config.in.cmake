@PACKAGE_INIT@

# Let CMake find included Find*.cmake files.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Required dependencies.
find_package(Eigen3 REQUIRED)
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)

# Optional dependencies.
set(RAYLIB_WITH_LAS @RAYLIB_WITH_LAS@)
if(RAYLIB_WITH_LAS)
  find_package(LASzip REQUIRED)
endif()
set(RAYLIB_WITH_QHULL @RAYLIB_WITH_QHULL@)
if(RAYLIB_WITH_QHULL)
  find_package(Qhull REQUIRED)
endif()
set(RAYLIB_WITH_TBB @RAYLIB_WITH_TBB@)
if(RAYLIB_WITH_TBB)
  find_package(TBB REQUIRED)
endif()
set(RAYLIB_WITH_TIFF @RAYLIB_WITH_TIFF@)
if(RAYLIB_WITH_TIFF)
  find_package(GeoTIFF REQUIRED)
  find_package(TIFF REQUIRED)
  pkg_check_modules(PROJ REQUIRED proj)
endif()

# Targets.
include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
check_required_components("@PROJECT_NAME@")
