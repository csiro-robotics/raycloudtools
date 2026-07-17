@PACKAGE_INIT@

# Required dependencies.
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# Optional dependencies.
set(RAYLIB_WITH_TIFF @RAYLIB_WITH_TIFF@)
if(RAYLIB_WITH_TIFF)
  find_package(PROJ REQUIRED)
endif()
set(RAYLIB_WITH_QHULL @RAYLIB_WITH_QHULL@)
if(RAYLIB_WITH_TIFF)
  find_package(Qhull REQUIRED)
endif()
set(RAYLIB_WITH_TBB @RAYLIB_WITH_TBB@)
if(RAYLIB_WITH_TIFF)
  find_package(TBB REQUIRED)
endif()

# Targets.
include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
check_required_components("@PROJECT_NAME@")
