if (!MSVC)
add_compile_options(-Wno-deprecated -Wno-write-strings -Wno-unused-result)
endif()
add_definitions(-DNDEBUG )
include_directories(../raylib/include)
find_package(libLAS REQUIRED)
find_package(libnabo REQUIRED)

# qhull
find_path(TRAVERSIBILITY_QHULL_INCLUDE_DIR libqhullcpp/Qhull.h)
if (NOT TRAVERSIBILITY_QHULL_INCLUDE_DIR)
   message(FATAL_ERROR "Unable to find libqhullcpp/Qhull.h")
endif()
find_library(TRAVERSIBILITY_LIBQHULLCPP qhullcpp libqhullcpp)
if (NOT TRAVERSIBILITY_LIBQHULLCPP)
   message(FATAL_ERROR "You need to have the qhull cpp library installed.")
endif()
SET(qhullDir /home/parallels/code/qhull-2015.2)
INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
INCLUDE_DIRECTORIES(${qhullDir}/src)
LINK_DIRECTORIES(${qhullDir}/build)
INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
INCLUDE_DIRECTORIES(include)
SET(qhullLibs qhullcpp qhull_r)


set(TARGETS
	raycreate
	raydenoise
	raysmooth
	raywrap
	rayconvert
	raydecimate
	raydiff
	raytranslate
	rayrotate
  raysplit
  raytransients
)
foreach(TARGET ${TARGETS})
	add_executable(${TARGET} ${TARGET}.cpp)
endforeach(TARGET)

foreach(TARGET ${TARGETS})
	target_link_libraries(${TARGET} raylib ${libLAS_LIBRARIES} ${libnabo_LIBRARIES} ${qhullLibs})
	set_target_properties(${TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../bin)
	set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${TARGET})
	install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
endforeach(TARGET)
