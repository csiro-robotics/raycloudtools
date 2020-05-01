# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Author: Kazys Stepanas

# This is a generalised find package helper script to be used when no <package>-config.cmake file is available. It
# ensures found libraries are always defined as import targets and resolves bindings for debug and release libraries.

include(CMakeParseArguments)

macro(_target_name IMPORT_TARGET)
  # Find include directories.
  string(TOUPPER TARGET_NAME "${IMPORT_TARGET}")
  string(REGEX REPLACE "[!@#\\$%\\^\\&\\*\\.\\(\\):;]" "_" TARGET_NAME "${IMPORT_TARGET}")
endmacro(_target_name)

###
function(ras_package_import_target IMPORT_TARGET)
  cmake_parse_arguments(PIT "INTERFACE;SHARED" "" "" ${ARGN})
  set(LIBRARY_ARGS UNKNOWN)
  if(PIT_INTERFACE)
    set(LIBRARY_ARGS INTERFACE)
  endif(PIT_INTERFACE)
  if(PIT_SHARED)
    set(LIBRARY_ARGS SHARED)
  endif(PIT_SHARED)

  add_library(${IMPORT_TARGET} ${LIBRARY_ARGS} IMPORTED)
endfunction(ras_package_import_target)

###
function(ras_package_import_include_dirs IMPORT_TARGET)
  _target_name(${IMPORT_TARGET})
  # Find include directories.
  find_path(${TARGET_NAME}_INCLUDE_DIR ${ARGN})

  if(PIM_DEBUG)
    message("+++++++++++++++++++++++++++")
    message("find_path(${TARGET_NAME}_INCLUDE_DIR ${ARGN})")
    message("${TARGET_NAME}_INCLUDE_DIR: ${${TARGET_NAME}_INCLUDE_DIR}")
  endif(PIM_DEBUG)

  if(${TARGET_NAME}_INCLUDE_DIR)
    set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${${TARGET_NAME}_INCLUDE_DIR}")
  endif(${TARGET_NAME}_INCLUDE_DIR)
endfunction(ras_package_import_include_dirs)


###
function(ras_package_import_library IMPORT_TARGET CONFIG)
  _target_name(${IMPORT_TARGET})
  string(TOUPPER "${CONFIG}" CONFIG_UPPER)

  find_library(${TARGET_NAME}_LIBRARY_${CONFIG_UPPER} ${ARGN})

  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
    find_library(${TARGET_NAME}_LIBRARY_DLL_${CONFIG_UPPER} ${IMPORT_TARGET})
  endif(WIN32)

  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".pdb")
    find_library(${TARGET_NAME}_LIBRARY_PDB_${CONFIG_UPPER} ${IMPORT_TARGET})
  endif(WIN32)

  if(${TARGET_NAME}_LIBRARY_${CONFIG_UPPER})
    set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${CONFIG_UPPER})
    # Check if the DLL is defined. If, so register IMPORTED_IMPLIB_<CONFIG> and IMPORTED_LOCATION_<CONFIG> otherwise,
    # just the later.
    if(${TARGET_NAME}_LIBRARY_DLL_${CONFIG_UPPER})
      set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY
        IMPORTED_IMPLIB_${CONFIG} "${${TARGET_NAME}_LIBRARY_${CONFIG_UPPER}}"
        IMPORTED_LOCATION_${CONFIG} "${${TARGET_NAME}_LIBRARY_DLL_${CONFIG_UPPER}}"
      )
    else(${TARGET_NAME}_LIBRARY_DLL_${CONFIG_UPPER})
      set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY
        IMPORTED_LOCATION_${CONFIG} "${${TARGET_NAME}_LIBRARY_${CONFIG_UPPER}}"
      )
    endif(${TARGET_NAME}_LIBRARY_DLL_${CONFIG_UPPER})
  endif(${TARGET_NAME}_LIBRARY_${CONFIG_UPPER})
endfunction(ras_package_import_library)

###
function(ras_package_import_include_directories IMPORT_TARGET)
  set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ARGN})
endfunction(ras_package_import_include_directories)

###
function(ras_package_import_link_libraries IMPORT_TARGET)
  set_property(TARGET ${IMPORT_TARGET} APPEND PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES ${ARGN})
endfunction(ras_package_import_link_libraries)

###
function(ras_package_import_validate IMPORT_TARGET)
  if(NOT IMPORT_TARGET)
    message(SEND_ERROR "${IMPORT_TARGET}: import target not defined")
    return()
  endif(NOT IMPORT_TARGET)

  get_target_property(TARGET_TYPE ${IMPORT_TARGET} TYPE)

  get_target_property(INCLUDE_DIRS ${IMPORT_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
    get_target_property(LIBRARY_BASE ${IMPORT_TARGET} IMPORTED_LOCATION)
    get_target_property(LIBRARY_RELEASE ${IMPORT_TARGET} IMPORTED_LOCATION_Release)
    get_target_property(LIBRARY_DEBUG ${IMPORT_TARGET} IMPORTED_LOCATION_Debug)

    # if(NOT LIBRARY_BASE AND NOT LIBRARY_RELEASE AND NOT LIBRARY_DEBUG)
    #   message(SEND_ERROR "${IMPORT_TARGET}: no import library found")
    # endif(NOT LIBRARY_BASE AND NOT LIBRARY_RELEASE AND NOT LIBRARY_DEBUG)
  endif(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
endfunction(ras_package_import_validate)

###
function(ras_package_import_define PACKAGE_NAME)
  cmake_parse_arguments(PID "" "" "TARGETS" ${ARGN})

  set(LOCAL_INCLUDE_DIRS)
  set(LOCAL_LIBRARIES)
  set(INCLUDE_NOT_FOUND FALSE)
  set(LIBRARY_NOT_FOUND FALSE)

  foreach(IMPORT_TARGET ${PID_TARGETS})
    _target_name(${IMPORT_TARGET})
    ras_package_import_validate(${IMPORT_TARGET})

    get_target_property(TARGET_TYPE ${IMPORT_TARGET} TYPE)

    get_target_property(INCLUDE_DIRS ${IMPORT_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    list(APPEND LOCAL_INCLUDE_DIRS ${INCLUDE_DIRS})

    if(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
      if(WIN32)
        get_target_property(LIBRARY_RELEASE ${IMPORT_TARGET} IMPORTED_IMPLIB_RELEASE)
      else(WIN32)
        get_target_property(LIBRARY_RELEASE ${IMPORT_TARGET} IMPORTED_LOCATION_RELEASE)
      endif(WIN32)

      if(WIN32)
        get_target_property(LIBRARY_DEBUG ${IMPORT_TARGET} IMPORTED_IMPLIB_DEBUG)
      else(WIN32)
        get_target_property(LIBRARY_DEBUG ${IMPORT_TARGET} IMPORTED_LOCATION_DEBUG)
      endif(WIN32)

      if(LIBRARY_DEBUG AND LIBRARY_RELEASE)
        list(APPEND LOCAL_LIBRARIES debug ${LIBRARY_DEBUG} optimized ${LIBRARY_RELEASE})
      elseif(LIBRARY_DEBUG)
        list(APPEND LOCAL_LIBRARIES ${LIBRARY_DEBUG})
      elseif(LIBRARY_RELEASE)
        list(APPEND LOCAL_LIBRARIES ${LIBRARY_RELEASE})
      else()
        set(APPEND LOCAL_LIBRARIES "${TARGET_NAME}-NOT_FOUND")
      endif()
    endif(NOT TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
  endforeach(IMPORT_TARGET)

  if(LOCAL_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES LOCAL_INCLUDE_DIRS)
  endif(LOCAL_INCLUDE_DIRS)
  if(LOCAL_LIBRARIES)
    list(REMOVE_DUPLICATES LOCAL_LIBRARIES)
  endif(LOCAL_LIBRARIES)

  set(${PACKAGE_NAME}_INCLUDE_DIRS "${LOCAL_INCLUDE_DIRS}" CACHE PATH "${PACKAGE_NAME} include directories" FORCE)
  set(${PACKAGE_NAME}_LIBRARIES "${LOCAL_LIBRARIES}" CACHE PATH "${PACKAGE_NAME} libraries" FORCE)
endfunction(ras_package_import_define)
