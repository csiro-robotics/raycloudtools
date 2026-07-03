@PACKAGE_INIT@
find_package(libnabo REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
check_required_components("@PROJECT_NAME@")
