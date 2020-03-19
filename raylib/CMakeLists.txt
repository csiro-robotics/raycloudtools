cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_SUPPRESS_REGENERATION true)
project("raylib")

#if (!MSVC)
	add_compile_options(-g) # -O3 -Wall)
#else()
#	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
#endif()

add_subdirectory(raylib/src)
