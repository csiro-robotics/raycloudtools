# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

include(CMakeParseArguments)

#-------------------------------------------------------------------------------
# ras_project(
#     [VERSION <major.minor.patch>]
#     [CXX_STD <14,17,...>])
#     [DEBUG_POSTFIX <postfix>]
#     [NAMESPACE <project-namespace>]
#     [PREFIX <include-prefix>]
#     [SRC_DIR <src-dir>]
# )
#
# Setup common project details and configuration. ras_project() is to be used after the standard CMAKE project()
# directive in order to configure standard aspects of a RAS project. In particular this covers the following
# aspects:
#
# **Version number**
# The version number is specified by the keyword variable "VERSION" followed by a version number in the format
# "<major>.<minor>.<patch>". This version number is split and extracted into the following variables:
# - PACKAGE_VERSION set to exactly match the passed in VERSION string.
# - PACKAGE_VERSION_MAJOR the major version part of the VERSION string.
# - PACKAGE_VERSION_MINOR the minor version part of the VERSION string.
# - PACKAGE_VERSION_PATCH the pathc version part of the VERSION string.
# - Aliases: each of the items above are have an alias. The alias is is named by replacing the "PACKAGE" part of the
#   variables above with ${CMAKE_PROJECT_NAME}. For example, if the project name is "raycloud", then the major version
#   number is available via "PACKAGE_VERSION_MAJOR" and "raycloudtools_VERSION_MAJOR".
#
# **C++ standard**
# The C++ standard is set by CXX_STD and can be one of [11, 14, 17]. The default is currently 14 and is used when
# CXX_STD is not present. The default standard and the list of valid standards will change as new C++ standards are
# more widely adopted.
#
# **Debug postfix**
# The CMAKE_DEBUG_POSTFIX is set according to the specified value of DEBUG_POSTFIX or the default string "d" when
# DEBUG_POSTFIX is not present. This is added as a suffix to the library name for all static, shared and module
# libraries. This is important for differentiating debug and release libraries, especially on the Windows platform.
#
# **Export namespace**
# The export namespace is set by the NAMESPACE value. This simply sets the value of the PACKAGE_NAMESPACE variable
# which is used by export() commands in ras_add_library(), ras_add_executable() directives.
#
# **Include prefix**
# The PREFIX is used when installing header files in ras_add_library() as the first level directly after the include
# directory. That is, header files are installed to "include/${PREFIX}".
#
# **Src directory**
# The SRC_DIR is used when marshalling build files in ras_add_library() as the first level directly after the 
# project folder directory. That is, build files are marshalled into "{PROJECT_BINARY_DIR}/{SRC_DIR}/{LIB_NAME}"/
# Note this is an optional prefix and the directories will collapse to "{PROJECT_BINARY_DIR}/{LIB_NAME} otherwise"
#
# **General configuration**
# Additionally, the following values are configured for the project.
#
# - PACKAGE_EXPORT_LOCATION lib/cmake/${CMAKE_PROJECT_NAME} - target for CMake package export scripts.
# - CMAKE_CXX_STANDARD_REQUIRED TRUE
# - CMAKE_POSITION_INDEPENDENT_CODE ON
# - CMAKE_EXPORT_COMPILE_COMMANDS ON
# - CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
# - GLOBAL PROPERTY USE_FOLDERS ON
#
# **Additional script import**
# - The following RAS CMake scripts are also imported:
# - RasCompilerSetup is included to configure base line warnings and compiler switches.
# - RasClangTidy
# - RasAddExecutable
# - RasAddLibrary
# - RasAddTest
# - RasUtil
#-------------------------------------------------------------------------------
macro(ras_project)
  cmake_parse_arguments(RAS "" "CXX_STD;DEBUG_POSTFIX;NAMESPACE;PREFIX;SRC_DIR;VERSION" "" ${ARGN})

  set(PACKAGE_VERSION "0.0.0")

  if(RAS_VERSION)
    set(PACKAGE_VERSION "${RAS_VERSION}")
  endif(RAS_VERSION)

  set(RAS_VERSION_REGX "([^\.]*)\.([^\.]*)\.([^\.]*)")
  string(REGEX REPLACE "${RAS_VERSION_REGX}" "\\1" PACKAGE_VERSION_MAJOR "${RAS_VERSION}")
  string(REGEX REPLACE "${RAS_VERSION_REGX}" "\\2" PACKAGE_VERSION_MINOR "${RAS_VERSION}")
  string(REGEX REPLACE "${RAS_VERSION_REGX}" "\\3" PACKAGE_VERSION_PATCH "${RAS_VERSION}")
  unset(RAS_VERSION_REGX)

  set(${CMAKE_PROJECT_NAME}_VERSION "${PACKAGE_VERSION}")
  set(${CMAKE_PROJECT_NAME}_VERSION_MAJOR "${PACKAGE_VERSION_MAJOR}")
  set(${CMAKE_PROJECT_NAME}_VERSION_MINOR "${PACKAGE_VERSION_MINOR}")
  set(${CMAKE_PROJECT_NAME}_VERSION_PATCH "${PACKAGE_VERSION_PATCH}")

  # Setup global variables used in configuring libraries and export
  set(PACKAGE_EXPORT_LOCATION lib/cmake/${CMAKE_PROJECT_NAME})
  if(RAS_NAMESPACE)
    set(PACKAGE_NAMESPACE ${RAS_NAMESPACE}::)
  endif(RAS_NAMESPACE)
  if(RAS_PREFIX)
    set(PACKAGE_INCLUDE_PREFIX ${RAS_PREFIX})
  endif(RAS_PREFIX)
  if(RAS_SRC_DIR)
    SET(SRC_DIR "/${RAS_SRC_DIR}")
  endif(RAS_SRC_DIR)

  # C++ standards setup.
  if(MSVC) # Should really look for VS2019+
    set(CMAKE_CXX_STANDARD 17)
  else(MSVC)
    set(CMAKE_CXX_STANDARD 14)
  endif(MSVC)
  if(RAS_CXX_STD)
    set(CMAKE_CXX_STANDARD ${RAS_CXX_STD})
  endif(RAS_CXX_STD)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
  # Ensure -fPIC is added.
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # Export compile_commands.json
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  # Ensure debug libraries are built with different named to release builds. This is to address issues such as MSVC
  # having different debug and release runtime libraries. For a well setup API, one which hides resource allocation and
  # ensures symmetrical deallcation occurs from the same allocator, this won't be a problem, but the consistency is
  # useful.
  set(CMAKE_DEBUG_POSTFIX "d")
  if(RAS_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX "${RAS_DEBUG_POSTFIX}")
  endif(RAS_DEBUG_POSTFIX)
  # Marshall all binaries to the same directory. This is expecially useful on Windows when trying to run exectuables from
  # this project with shared libraries. Otherwise those shared libraries aren't on the path. Note that other package
  # binaries should be on the path already.
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

  # Allow the use of folders to group targets in supporting environments.
  # For example Visual Studio solution folders.
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Manage compiler.
  # Use CMAKE_MODULE_PATH and include(compilerSetup) if compilerSetup.cmake is moved.
  include(RasCompilerSetup)
  include(RasClangTidy)
  include(RasUtil)
  include(RasAddLibrary)
  include(RasAddExecutable)
  include(RasAddTest)
endmacro(ras_project)
