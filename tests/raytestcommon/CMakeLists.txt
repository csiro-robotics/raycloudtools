# Setup of GTEST changed at CMake 3.5.
cmake_minimum_required(VERSION 3.5)

# Eigen required to support some tests
find_package(Eigen3 QUIET)

set(SOURCES
  raytestutils.cpp
  raytestutils.h
)

add_library(raytestcommon ${SOURCES})

set_target_properties(raytestcommon PROPERTIES FOLDER tests)

target_include_directories(raytestcommon
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests>
)

target_include_directories(raytestcommon SYSTEM
  PRIVATE
)

target_link_libraries(raytestcommon PUBLIC raytools ray rayutil)

if(Eigen3_FOUND)
  target_link_libraries(raytestcommon PUBLIC Eigen3::Eigen)
endif(Eigen3_FOUND)

source_group("source" REGULAR_EXPRESSION ".*$")
