# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230

set(SOURCES
  environment.cpp
  raytests.cpp
  testmain.cpp
)
add_executable(raytest ${SOURCES})
target_include_directories(raytest
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/raylib>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
)
target_include_directories(raytest SYSTEM
  PRIVATE
    "${GTEST_INCLUDE_DIRS}"
)
target_link_libraries(raytest
  PUBLIC
    raylib
    GTest::Main
)
target_link_libraries(raytest
  PRIVATE
    Eigen3::Eigen
)
gtest_discover_tests(raytest
  EXTRA_ARGS ${CMAKE_BINARY_DIR}
)
