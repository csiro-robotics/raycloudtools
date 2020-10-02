# Setup of GTEST changed at CMake 3.5.
cmake_minimum_required(VERSION 3.5)

find_package(GTest)
# Eigen required to support some tests - NDT in particular
find_package(Eigen3 QUIET)

set(SOURCES
  raytests.cpp
  testmain.cpp
)

add_executable(raytest ${SOURCES})
set_target_properties(raytest PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(raytest PROPERTIES FOLDER tests)

target_include_directories(raytest
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/raylib>    
)

target_include_directories(raytest SYSTEM
  PRIVATE
    "${GTEST_INCLUDE_DIRS}"
)

target_link_libraries(raytest PUBLIC raylib GTest::Main)

if(Eigen3_FOUND)
  target_link_libraries(raytest PRIVATE Eigen3::Eigen)
endif(Eigen3_FOUND)

add_test(NAME raytest COMMAND raytest --gtest_output=xml:test-reports/)

source_group("source" REGULAR_EXPRESSION ".*$")
