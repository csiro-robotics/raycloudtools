# Copyright (c) 2019-2020
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Authors: Gavin Catt, Kazys Stepanas

include(GenerateExportHeader)
include(RasUtil)
include(RasSourceGroupExtensions)

# Helper macro to add a library for adding a library in a standard way
# Usage:
# ras_add_library(<lib_name>
#   [PRIVATE_LIBRARY]
#   [TYPE <INTERFACE|STATIC|SHARED|MODULE>]
#   [INCLUDE [PUBLIC|PRIVATE|PUBLIC_SYSTEM|PRIVATE_SYSTEM]] <header1> ...
#   [LIBS [PUBLIC|PRIVATE|INTERFACE]] <library1> ...
#   [INCLUDE_PREFIX <sub_folder>]
#   [PROJECT_FOLDER <project_folder>]
#   [PUBLIC_HEADERS <header1> <header2> ...]
#   SOURCES <source1> <source2> ...
#   [GENERATED [PUBLIC|PRIVATE] <gen1> <gen2> ...]
#   [PCH <pch-header>]
# )
#
# TYPE may be one of INTERFACE, STATIC, SHARED or MODULE and defines the library type. All except INTERFACE represent
# standard CMake library types. INTERFACE identifies the library as an interface library; one which contains no
# compilable symbols. The library is build as a static library and dummy symbols are added to avoid compiler warnings.
# Ommitting the TYPE has the same effect as add_library() without specifying the type; i.e., default to STATIC unless
# BUILD_SHARED is true, in which case the type is SHARED.
#
# The INCLUDE and LIBS arguments support modifiers PUBLIC|PRIVATE with the addition of PUBLIC_SYSTEM and PRIVATE_SYSTEM
# modifier for INCLUDE directories. These control how these are passed to the target_include_directories() and
# target_link_libraries() commands. The PUBLIC modifier is implied and may be omitted. Both the PUBLIC and PRIVATE
# modifiers are essentially modal and affect everything that follows. The two SYSTEM  modifiers ensure include
# directories are added as SYSTEM includes and avoid compiler warnings from those external header files.
#
# Note that each of these modifier combinations can only appear once for INCLUDE and once for LIBS. Some examples below
# help illustrate this point.
#
# INCLUDE_PREFIX is used to identify the installation sub-folder for include directories. See PACKAGE_INCLUDE_PREFIX
# below.
#
# PROJECT_FOLDER is used to group libraries within IDEs which support sub-project structures such as Visual Studio
# (full version) or Xcode. This is ignored for most build systems.
#
# PUBLIC_HEADERS idnetifies the list of public header files for this library. These are included in the library
# installation as part of the public API.
#
# SOURCES identifies the list of source files for the library. Note: generated source file and header artefacts should
# not be included in this SOURCES list, but should be added using the GENERATED keyword argument.
#
# GENERATED identifies the list of source and header file artefacts which have been generated. These are expected to
# appear within the build directory and not within the source tree. GENERATED items support PUBLIC|PRIVATE modifiers
# in a similar fashion to INCLUDE and LIBS lists, however the implied modifier is PRIVATE. PUBLIC generated items are
# marshalled as part of the installation process, while PRIVATE items remain hidden.
#
# Examples:
# # Allowed usage.
# ras_add_library(mylib
#     INCLUDE PUBLIC <public1> <public2> ... PUBLIC_SYSTEM <public_system1> <public_system2> ...
#             PRIVATE <private1> <private2> PUBLIC_SYSTEM <private_system1> <private_system2>
# )
#
# Invalid usage (PUBLIC used twice and PUBLIC SYSTEM)
# ras_add_library(mylib
#     INCLUDE PUBLIC <public1> <public2> ... SYSTEM <public_system1> <public_system2> ...
#             PRIVATE <private1> <private2> SYSTEM <private_system1> <private_system2>
#             PUBLIC <public3> <public4> ... SYSTEM <public_system3> <public_system4> ... # Not allowed
# )
#
# Additionally, some global variables affect the library installation.
# PACKAGE_EXPORT_LOCATION - A project wide variable which defines where to install the package configuration files to.
#
# PACKAGE_NAMESPACE - defines the namespace prefix appended to the library name when exported to for use in by other
# projects. For example, the library 'core' library is exported as 'ras::core' with the PACKAGE_NAMESPACE 'ras'.
# May be empty
#
# PACKAGE_INCLUDE_PREFIX - a path appended to the include installation directory. Public headers are exported to the
# path: <install-dir>/include[/PACKAGE_INCLUDE_PREFIX]/<INCLUDE_PREFIX-argument>.
#
# PRIVATE_LIBRARY - indicates the library is private for build time only and should not be installed.
function(ras_add_library LIB_NAME)
  # Multi-value arguments to parse. Too long for a single string.
  cmake_parse_arguments(ARG "PRIVATE_LIBRARY" "INCLUDE_PREFIX;PCH;PROJECT_FOLDER;TYPE"
    "GENERATED;INCLUDE;LIBS;PUBLIC_HEADERS;SOURCES"
    ${ARGN})

  # Resolve library type.
  if(ARG_TYPE)
    # Validate library type.
    set(LIB_TYPES INTERFACE MODULE SHARED STATIC)
    list(FIND LIB_TYPES ${ARG_TYPE} LIB_TYPE_INDEX)
    if(LIB_TYPE_INDEX LESS 0)
      message(SEND_ERROR "Unknown library type ${ARG_TYPE}. Must be one of ${LIB_TYPES}")
    endif(LIB_TYPE_INDEX LESS 0)
  endif(ARG_TYPE)

  # Split includes into unspecified, public, public system, private, private system.
  ras_split_list(ARG_INCLUDE "PUBLIC;PRIVATE;PUBLIC_SYSTEM;PRIVATE_SYSTEM" ${ARG_INCLUDE})
  # Merge ARG_INCLUDE and ARG_INCLUDE_PUBLIC
  list(APPEND ARG_INCLUDE_PUBLIC ${ARG_INCLUDE})

  if(AWL_DEBUG)
    message("ARG_INCLUDE: ${ARG_INCLUDE}")
    message("ARG_INCLUDE_PUBLIC: ${ARG_INCLUDE_PUBLIC}")
    message("ARG_INCLUDE_PRIVATE: ${ARG_INCLUDE_PRIVATE}")
    message("ARG_INCLUDE_PUBLIC_SYSTEM: ${ARG_INCLUDE_PUBLIC_SYSTEM}")
    message("ARG_INCLUDE_PRIVATE: ${ARG_INCLUDE_PRIVATE}")
  endif(AWL_DEBUG)

  if(AWL_DEBUG)
    message("ARG_LIBS(in): ${ARG_LIBS}")
  endif(AWL_DEBUG)
  # # Split libs into unspecified, public, private.
  ras_split_list(ARG_LIBS "PUBLIC;PRIVATE;INTERFACE" ${ARG_LIBS})
  # Merge ARG_LIBS and ARG_LIBS_PUBLIC
  list(APPEND ARG_LIBS_PUBLIC ${ARG_LIBS})

  if(AWL_DEBUG)
    message("ARG_LIBS: ${ARG_LIBS}")
    message("ARG_LIBS_PUBLIC: ${ARG_LIBS_PUBLIC}")
    message("ARG_LIBS_PRIVATE: ${ARG_LIBS_PRIVATE}")
    message("ARG_LIBS_INTERFACE: ${ARG_LIBS_INTERFACE}")
  endif(AWL_DEBUG)

  if(AWL_DEBUG)
    message("ARG_GENERATED(in): ${ARG_GENERATED}")
  endif(AWL_DEBUG)
  # Split GENERATED items into unspecified, public, private.
  ras_split_list(ARG_GENERATED "PUBLIC;PRIVATE" ${ARG_GENERATED})
  # Merge ARG_GENERATED and ARG_GENERATED_PRIVATE
  list(APPEND ARG_GENERATED_PRIVATE ${ARG_GENERATED})

  # Recombine ARG_GENERATED as the collation of ARG_GENERATED_PUBLIC and ARG_GENERATED_PRIVATE
  set(ARG_GENERATED ${ARG_GENERATED_PUBLIC} ${ARG_GENERATED_PRIVATE})

  if(AWL_DEBUG)
    message("ARG_GENERATED: ${ARG_GENERATED}")
    message("ARG_GENERATED_PUBLIC: ${ARG_GENERATED_PUBLIC}")
    message("ARG_GENERATED_PRIVATE: ${ARG_GENERATED_PRIVATE}")
  endif(AWL_DEBUG)

  # In the event of there being no module sources (.cpp, etc.)
  # Then we need to build an INTERFACE library
  # set ARG_INTERFACE if interface lib

  # Catkin does NOT support INTERFACE libs as it attempts to link the targets as straight up libs
  # This means we need to create an actual library, to do this we just generate a new lib with some dummy symbols
  if(ARG_TYPE STREQUAL INTERFACE)
    message(STATUS "Library '${LIB_NAME}' is an interface library; adding dummy symbols.")
    set(ARG_TYPE STATIC)
    # Generate dummy files. Only do so if they don't already exist.
    set(DUMMY_LIB_NAME ${LIB_NAME})
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Dummy.h")
      configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/Dummy.h.in"
        Dummy.h
      )
    endif(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Dummy.h")
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Dummy.cpp")
      configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/Dummy.cpp.in"
        Dummy.cpp
      )
    endif(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Dummy.cpp")
      # Add them to the build
    list(APPEND ARG_SOURCES "Dummy.cpp" "Dummy.h")
  endif(ARG_TYPE STREQUAL INTERFACE)

  # Add the lib
  add_library(${LIB_NAME} ${ARG_TYPE} ${ARG_PUBLIC_HEADERS} ${ARG_SOURCES} ${ARG_GENERATED})
  # Export the lib to our target list
  set(PROJECT_LIBRARY_TARGETS "${PROJECT_LIBRARY_TARGETS};${LIB_NAME}" CACHE INTERNAL PROJECT_LIBRARY_TARGETS)

  # Generate an export header for DLL usage.
  generate_export_header(${LIB_NAME})
  list(APPEND ARG_GENERATED_PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}_export.h")

  # Setup solution display for the target.
  # Set the solution folder for the target (optional).
  if(ARG_PROJECT_FOLDER)
    set_target_properties(${LIB_NAME} PROPERTIES FOLDER ${ARG_PROJECT_FOLDER})
  endif(ARG_PROJECT_FOLDER)

  # Identify public headers. Will be installed. Generated headers to be installed should be included in the list.
  # set_target_properties(${LIB_NAME} PROPERTIES PUBLIC_HEADER "${ARG_GENERATED_PUBLIC};${ARG_PUBLIC_HEADERS}")

  # Setup target versioning.
  if(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)
    set_target_properties(${LIB_NAME} PROPERTIES
      VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
    )
  endif(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)

  # Add include directories to the target
  target_include_directories(${LIB_NAME}
    PUBLIC
      # Public include directories.
      # Include either PROJECT_BINARY_DIR or CMAKE_CURRENT_BINARY_DIR during build to access the config and export
      # headers. The former supports using relative include paths such as "#include <ras/core/Config.h>"", while the
      # latter is for an unprefixed approach; "#include <RasConfig.h>"
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
      # # Also include the directory up for external targets correctly including the generated headers.
      # $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
      # Include the project directory to get access to other libraries in the project.
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
      # Include the installation directory for import.
      $<INSTALL_INTERFACE:include>
  )
    
  if(AWL_DEBUG)
    message("INCLUDE:")
    message("PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")
    message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")  
  endif(AWL_DEBUG)

  # Includes
  target_include_directories(${LIB_NAME}
    PUBLIC ${ARG_INCLUDE_PUBLIC}
    PRIVATE ${ARG_INCLUDE_PRIVATE}
  )
  target_include_directories(${LIB_NAME} SYSTEM
    PUBLIC ${ARG_INCLUDE_PUBLIC_SYSTEM}
    PRIVATE ${ARG_INCLUDE_PRIVATE_SYSTEM}
  )

  list(FIND ARG_LIBS_PUBLIC SYSTEM SEARCH_INDEX)
  if(SEARCH_INDEX GREATER -1)
    message(FATAL_ERROR "Adding PUBLIC SYSTEM link for ${LIB_NAME}. SYSTEM should only be used for includes.")
  endif(SEARCH_INDEX GREATER -1)

  list(FIND ARG_LIBS_PRIVATE SYSTEM SEARCH_INDEX)
  if(SEARCH_INDEX GREATER -1)
    message(FATAL_ERROR "Adding PRIVATE SYSTEM link for ${LIB_NAME}. SYSTEM should only be used for includes.")
  endif(SEARCH_INDEX GREATER -1)

  # Link deps
  target_link_libraries(${LIB_NAME}
    PUBLIC ${ARG_LIBS_PUBLIC} ${CMAKE_DL_LIBS}
    PRIVATE ${ARG_LIBS_PRIVATE}
    INTERFACE ${ARG_LIBS_INTERFACE}
  )

  # Enable clang-tidy
  ras_clang_tidy_target(${LIB_NAME} EXCLUDE_MATCHES ".*\\.in($|\\..*)")

  if(NOT ARG_PRIVATE_LIBRARY)
    # Setup installation.
    set(INCLUDE_INSTALL_DIR "include/${PACKAGE_INCLUDE_PREFIX}/${ARG_INCLUDE_PREFIX}")

    # Header installation
    # Note we want to preserve there folder structure from the root so we can't use PUBLIC_HEADER DESTINATION
    foreach (file ${ARG_PUBLIC_HEADERS})
      get_filename_component(dir ${file} DIRECTORY)
      install(FILES ${file} DESTINATION "${INCLUDE_INSTALL_DIR}/${dir}")
    endforeach(file)

    foreach (file ${ARG_GENERATED_PUBLIC})
      # Get the generated file directory relative to CMAKE_CURRENT_BINARY_DIR
      get_filename_component(dir "${file}" DIRECTORY)
      get_relative_path(dir "${dir}" "${CMAKE_CURRENT_BINARY_DIR}")
      install(FILES ${file} DESTINATION "${INCLUDE_INSTALL_DIR}/${dir}")
    endforeach(file)

    # Binary installation
    install(TARGETS ${LIB_NAME} EXPORT ${CMAKE_PROJECT_NAME}-targets
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      RUNTIME DESTINATION bin
      # Defines the base include path
      INCLUDES DESTINATION "include"
    )

    # Install PDB files (MSVC)
    if(MSVC)
      if(BUILD_SHARED_LIBS OR ARG_TYPE STREQUAL "MODULE" OR ARG_TYPE STREQUAL "SHARED")
        install(FILES $<TARGET_PDB_FILE:${LIB_NAME}> DESTINATION bin OPTIONAL)
      endif(BUILD_SHARED_LIBS OR ARG_TYPE STREQUAL "MODULE" OR ARG_TYPE STREQUAL "SHARED")
    endif(MSVC)

    export(EXPORT ${CMAKE_PROJECT_NAME}-targets
      FILE "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}-targets.cmake"
      NAMESPACE ${PACKAGE_NAMESPACE}
    )
  endif(NOT ARG_PRIVATE_LIBRARY)

  # Setup folder display with the target for Visual Studio. This should always be done to match
  # the on disk layout of the source files.
  ras_make_source_groups(GENERATED ${ARG_GENERATED} SOURCE ${ARG_SOURCES} ${ARG_PUBLIC_HEADERS})
endfunction(ras_add_library)
