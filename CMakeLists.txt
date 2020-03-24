cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_SUPPRESS_REGENERATION true)
project("raycloudtools")

add_subdirectory(raylib/src)
add_subdirectory(src)
