if (!MSVC)
	add_compile_options(-Wno-deprecated -Wno-write-strings -Wno-unused-result -std=c++11)
endif()
add_definitions(-DNDEBUG -DUNORDERED -DHAVE_UNORDERED_MAP)

file(GLOB_RECURSE SOURCE_FILES
	${CMAKE_SOURCE_DIR}/raylib/src/*.cpp
)

file(GLOB_RECURSE INCLUDE_FILES
	${CMAKE_SOURCE_DIR}/raylib/include/*.h
)

add_compile_options(-std=c++11 -Wall -g)
find_package(Eigen3 REQUIRED)
find_package(libLAS REQUIRED)

add_library(raylib ${SOURCE_FILES})
set_property(TARGET raylib PROPERTY POSITION_INDEPENDENT_CODE ON)

if (BUILD_SHARED_LIBS)
	target_compile_definitions(raylib PRIVATE "COMPILE_AS_DLL")
endif()

target_include_directories(raylib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/raylib/include>
  $<INSTALL_INTERFACE:include/raylib>
  ${EIGEN3_INCLUDE_DIR}
)
set_target_properties(raylib PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
)

if(MSVC)
   target_compile_options(raylib PRIVATE /Zi)

   # Tell linker to include symbol data
    set_target_properties(raylib PROPERTIES 
        LINK_FLAGS "/INCREMENTAL:NO /DEBUG /OPT:REF /OPT:ICF"
    )

    # Set file name & location
    set_target_properties(raylib PROPERTIES 
        COMPILE_PDB_NAME raylib 
        COMPILE_PDB_OUTPUT_DIR ${CMAKE_BINARY_DIR}
    )
endif()

install(FILES ${INCLUDE_FILES} DESTINATION include/raylib)

if (MSVC)
	foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
		install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib/${OUTPUTCONFIG} DESTINATION lib/raylib)
	endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
else()
	install(TARGETS raylib EXPORT raylib-targets
		ARCHIVE DESTINATION lib/raylib
		LIBRARY DESTINATION lib/raylib
		RUNTIME DESTINATION lib/raylib)
	install(EXPORT raylib-targets DESTINATION lib/cmake/raylib)
	install(FILES ${CMAKE_SOURCE_DIR}/raylib/src/raylib-config.cmake DESTINATION lib/cmake/raylib)
endif(MSVC)
