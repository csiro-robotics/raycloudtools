include(RasPackageImport)

# Import Qhull Cpp target
ras_package_import_target(Qhull::qhullcpp)
ras_package_import_include_dirs(Qhull::qhullcpp Qhull.h PATH_SUFFIXES libqhullcpp)
ras_package_import_library(Qhull::qhullcpp Debug qhullcpp libqhullcpp)
ras_package_import_library(Qhull::qhullcpp Release qhullcpp libqhullcpp)
ras_package_import_define(Qhull TARGETS Qhull::qhullcpp)

# Import Qhull static reentrant target
ras_package_import_target(Qhull::qhullstatic_r)
ras_package_import_include_dirs(Qhull::qhullstatic_r libqhull_r.h PATH_SUFFIXES libqhull_r)
ras_package_import_library(Qhull::qhullstatic_r Debug qhull_r libqhull_r)
ras_package_import_library(Qhull::qhullstatic_r Release qhull_r libqhull_r)
ras_package_import_define(Qhull TARGETS Qhull::qhullstatic_r)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull REQUIRED_VARS
  Qhull_INCLUDE_DIRS)
mark_as_advanced(Qhull_INCLUDE_DIRS)