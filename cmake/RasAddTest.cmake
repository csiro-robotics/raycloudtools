# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Authors: Gavin Catt, Kazys Stepanas

include(RasUtil)
include(RasSourceGroupExtensions)

# Helper macro to add a test for ras to the build tree with some pre-defined inputs
# Usage:
# ras_add_test(<test_name>
#   [INCLUDE [PUBLIC|PRIVATE|PUBLIC_SYSTEM|PRIVATE_SYSTEM]] <header1> ...
#   [LIBS [PUBLIC|PRIVATE]] <library1> ...
#   [PROJECT_FOLDER] <sub_folder>
#   SOURCES <source1> <source2> ...
#   [GENERATED [PUBLIC|PRIVATE] <gen1> <gen2> ...]
# )
#
# The behaviour of various keyword arguments is identical to that used for ras_add_library. There are simply fewer
# such keyword arguments due to more limited requirements for executable targets.
function(ras_add_test TEST_NAME)
  # Multi-value arguments to parse. Too long for a single string.
  cmake_parse_arguments(ARG "" "PROJECT_FOLDER" "INCLUDE;LIBS;SOURCES" ${ARGN})

  # Split includes into unspecified, public, public system, private, private system.
  ras_split_list(ARG_INCLUDE "PUBLIC;PRIVATE;PUBLIC_SYSTEM;PRIVATE_SYSTEM" ${ARG_INCLUDE})
  # Merge ARG_INCLUDE and ARG_INCLUDE_PUBLIC
  list(APPEND ARG_INCLUDE_PUBLIC ${ARG_INCLUDE})

  # # Split libs into unspecified, public, private.
  ras_split_list(ARG_LIBS "PUBLIC;PRIVATE" ${ARG_LIBS})
  # Merge ARG_LIBS and ARG_LIBS_PUBLIC
  list(APPEND ARG_LIBS_PUBLIC ${ARG_LIBS})

  # Split GENERATED items into unspecified, public, private.
  ras_split_list(ARG_GENERATED "PUBLIC;PRIVATE" ${ARG_GENERATED})
  # Merge ARG_GENERATED and ARG_GENERATED_PRIVATE
  list(APPEND ARG_GENERATED_PRIVATE ${ARG_GENERATED})

  # TODO(KS): make sure this -rdynamic addition works then move into the root CMakeLists file.
  # Add -rdyanmic
  set(CMAKE_ENABLE_EXPORTS ON)
  # set(CMAKE_C_FLAGS "-rdynamic")
  # set(CMAKE_CXX_FLAGS "-rdynamic")
  # set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")

  add_executable(${TEST_NAME} ${ARG_SOURCES})

  # CMake does not automatically propagate CMAKE_DEBUG_POSTFIX to executables. We do so to avoid confusing link issues
  # which can would when building release and debug exectuables to the same path.
  set_target_properties(${TEST_NAME} PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

  # Setup solution display for the target.
  # Set the solution folder for the target (optional).
  if(ARG_FOLDER)
    set_target_properties(${TEST_NAME} PROPERTIES FOLDER ${ARG_FOLDER})
  endif(ARG_FOLDER)

  # # Setup target versioning.
  # if(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)
  #   set_target_properties(${TEST_NAME} PROPERTIES
  #     VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
  #   )
  # endif(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)

  # Add include directories to the target
  target_include_directories(${TEST_NAME}
    PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${PROJECT_SOURCE_DIR}"
  )
  target_include_directories(${TEST_NAME}
    PUBLIC ${ARG_INCLUDE_PUBLIC}
    PRIVATE ${ARG_INCLUDE_PRIVATE}
  )
  target_include_directories(${TEST_NAME} SYSTEM
    PUBLIC ${ARG_INCLUDE_PUBLIC_SYSTEM}
    PRIVATE ${ARG_INCLUDE_PRIVATE_SYSTEM}
  )

  # Link dependencies.
  target_link_libraries(${TEST_NAME}
    PUBLIC ${ARG_LIBS_PUBLIC}
    PRIVATE ${ARG_LIBS_PRIVATE} ${CMAKE_DL_LIBS}
    INTERFACE ${ARG_LIBS_INTERFACE}
  )

  # Define the tests to run.
  # With Google Test, we can either simply use the executable above as the command, or we can use --gtest_filter=XXX
  # to run individual or groups of tests (remember --gtest_filter supports wild-card matching).
  # Use --gtest_output=xml: to generate output which a CI server can interpret and report.
  add_test(NAME ${TEST_NAME} # Does not have to match the target name.
           COMMAND ${TEST_NAME} --gtest_output=xml:test-reports/)

  # Install PDB files (MSVC)
  if(MSVC)
    install(FILES $<TARGET_PDB_FILE:${TEST_NAME}> DESTINATION bin OPTIONAL)
  endif(MSVC)

  # Setup folder display with the target for Visual Studio. This should always be done to match
  # the on disk layout of the source files.
  ras_make_source_groups(GENERATED ${ARG_GENERATED} SOURCE ${ARG_SOURCES} ${ARG_PUBLIC_HEADERS})
endfunction(ras_add_test)
