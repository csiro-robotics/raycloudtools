# Copyright (c) 2019
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Authors: Gavin Catt, Kazys Stepanas

include(RasUtil)
include(RasSourceGroupExtensions)

# Helper function to add an application for ras to the build tree with some pre-defined inputs
# Usage:
# ras_add_executable(<app_name>
#   [INCLUDE [PUBLIC|PRIVATE|PUBLIC_SYSTEM|PRIVATE_SYSTEM]] <header1> ...
#   [LIBS [PUBLIC|PRIVATE|INTERFACE]] <library1> ...
#   [INCLUDE_PREFIX] <sub_folder>
#   [PROJECT_FOLDER] <sub_folder>
#   SOURCES <source1> <source2> ...
#   [GENERATED [PUBLIC|PRIVATE] <gen1> <gen2> ...]
# )
#
# The behaviour of various keyword arguments is identical to that used for ras_add_library. There are simply fewer
# such keyword arguments due to more limited requirements for executable targets.
function(ras_add_executable APP_NAME)
  # Multi-value arguments to parse. Too long for a single string.
  cmake_parse_arguments(ARG "" "INCLUDE_PREFIX;PROJECT_FOLDER" "GENERATED;INCLUDE;LIBS;SOURCES" ${ARGN})

  # Split includes into unspecified, public, public system, private, private system.
  ras_split_list(ARG_INCLUDE "PUBLIC;PRIVATE;PUBLIC_SYSTEM;PRIVATE_SYSTEM" ${ARG_INCLUDE})
  # Merge ARG_INCLUDE and ARG_INCLUDE_PUBLIC
  list(APPEND ARG_INCLUDE_PUBLIC ${ARG_INCLUDE})

  # # Split libs into unspecified, public, private.
  ras_split_list(ARG_LIBS "PUBLIC;PRIVATE" ${ARG_LIBS})
  # Merge ARG_LIBS and ARG_LIBS_PUBLIC
  list(APPEND ARG_LIBS_PUBLIC ${ARG_LIBS})

  # Split GENERATED items into unspecified, public, private.
  ras_split_list(ARG_GENERATED "PUBLIC;PRIVATE;INTERFACE" ${ARG_GENERATED})
  # Merge ARG_GENERATED and ARG_GENERATED_PRIVATE
  list(APPEND ARG_GENERATED_PRIVATE ${ARG_GENERATED})

  # Recombine ARG_GENERATED as the collation of ARG_GENERATED_PUBLIC and ARG_GENERATED_PRIVATE
  set(ARG_GENERATED ${ARG_GENERATED_PUBLIC} ${ARG_GENERATED_PRIVATE} ${ARG_GENERATED_INTERFACE})

  # TODO(KS): make sure this -rdynamic addition works then move into the root CMakeLists file.
  # Add -rdyanmic
  set(CMAKE_ENABLE_EXPORTS ON)
  # set(CMAKE_C_FLAGS "-rdynamic")
  # set(CMAKE_CXX_FLAGS "-rdynamic")
  # set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")

  add_executable(${APP_NAME} ${ARG_SOURCES} ${ARG_GENERATED})

  # CMake does not automatically propagate CMAKE_DEBUG_POSTFIX to executables. We do so to avoid confusing link issues
  # which can would when building release and debug exectuables to the same path.
  set_target_properties(${APP_NAME} PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

  # Setup solution display for the target.
  # Set the solution folder for the target (optional).
  if(ARG_PROJECT_FOLDER)
    set_target_properties(${APP_NAME} PROPERTIES FOLDER ${ARG_PROJECT_FOLDER})
  endif(ARG_PROJECT_FOLDER)

  # Setup target versioning.
  if(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)
    set_target_properties(${APP_NAME} PROPERTIES
      VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
    )
  endif(DEFINED ${CMAKE_PROJECT_NAME}_VERSION)

  # Add include directories to the target
  target_include_directories(${APP_NAME}
    PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${PROJECT_SOURCE_DIR}"
  )

  target_include_directories(${APP_NAME}
    PUBLIC ${ARG_INCLUDE_PUBLIC}
    PRIVATE ${ARG_INCLUDE_PRIVATE}
  )
  target_include_directories(${APP_NAME} SYSTEM
    PUBLIC ${ARG_INCLUDE_PUBLIC_SYSTEM}
    PRIVATE ${ARG_INCLUDE_PRIVATE_SYSTEM}
  )

  # Link dependencies.
  target_link_libraries(${APP_NAME}
    PUBLIC ${ARG_LIBS_PUBLIC}
    PRIVATE ${ARG_LIBS_PRIVATE} ${CMAKE_DL_LIBS}
    INTERFACE ${ARG_LIBS_INTERFACE}
  )

  # Enable clang-tidy
  ras_clang_tidy_target(${APP_NAME} EXCLUDE_MATCHES ".*\\.in($|\\..*)")

  # # Enable LeakTrack if required
  # ras_leak_track_target_enable(${APP_NAME} CONDITION RAYCLOUD_LEAK_TRACK)
  # # Suppress known false-positives/ignoreable leaks in third party libraries
  # ras_leak_track_suppress(${APP_NAME} CONDITION RAYCLOUD_LEAK_TRACK
  #   # /lib/x86_64-linux-gnu/libnss_db.so.2 (2 x 48 Bytes)
  #   "libnss_db"
  #   # /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (1 x 43 Bytes)
  #   "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >::_M_assign(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)"
  #   # /opt/ros/melodic/lib/libclass_loader.so (1 x 48 Bytes)
  #   "class_loader::impl::loadLibrary(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, class_loader::ClassLoader*)"
  # )

  set(INCLUDE_INSTALL_DIR include)
  if(PACKAGE_INCLUDE_PREFIX)
    set(INCLUDE_INSTALL_DIR "${INCLUDE_INSTALL_DIR}/${PACKAGE_INCLUDE_PREFIX}")
  endif(PACKAGE_INCLUDE_PREFIX)
  if(ARG_INCLUDE_PREFIX)
    set(INCLUDE_INSTALL_DIR "${INCLUDE_INSTALL_DIR}/${ARG_INCLUDE_PREFIX}")
  endif(ARG_INCLUDE_PREFIX)

  #Install
  install(TARGETS ${APP_NAME}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION "${INCLUDE_INSTALL_DIR}"
  )

  # Install PDB files (MSVC)
  if(MSVC)
    install(FILES $<TARGET_PDB_FILE:${APP_NAME}> DESTINATION bin OPTIONAL)
  endif(MSVC)

  # Setup folder display with the target for Visual Studio. This should always be done to match
  # the on disk layout of the source files.
  ras_make_source_groups(GENERATED ${ARG_GENERATED} SOURCE ${ARG_SOURCES} ${ARG_PUBLIC_HEADERS})
endfunction(ras_add_executable)
