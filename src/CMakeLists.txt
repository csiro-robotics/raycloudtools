if (!MSVC)
add_compile_options(-Wno-deprecated -Wno-write-strings -Wno-unused-result)
endif()
add_definitions(-DNDEBUG )
include_directories(../raylib/include)
find_package(libLAS REQUIRED)

set(TARGETS
	raycreate
	raydenoise
	raysmooth
	raywrap
	rayconvert
	raydecimate
	raydiff
	raytransform
)
foreach(TARGET ${TARGETS})
	add_executable(${TARGET} ${TARGET}.cpp)
endforeach(TARGET)

foreach(TARGET ${TARGETS})
	target_link_libraries(${TARGET} raylib ${libLAS_LIBRARIES})
	set_target_properties(${TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../bin)
	set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${TARGET})
	install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
endforeach(TARGET)
