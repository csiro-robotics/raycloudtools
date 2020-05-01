@PACKAGE_INIT@
# set(@PROJECT_NAME@_DIR "@PACKAGE_SOME_INSTALL_DIR@")

include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
check_required_components("@PROJECT_NAME@")

# Add public cmake modules.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
