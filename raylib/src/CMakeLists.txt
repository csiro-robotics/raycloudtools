
add_compile_options(-g -O3 -Wall -fopenmp)
add_definitions(-DNDEBUG -DUNORDERED -DHAVE_UNORDERED_MAP)

IF(USE_LAS)
  message("raylib: USE_LAS enabled")
add_definitions(-DUSE_LAS)
ENDIF(USE_LAS)

IF(USE_QHULL)
  message("raylib: USE_QHULL enabled")
  add_definitions(-DUSE_QHULL)
ENDIF(USE_QHULL)

IF(USE_ROS)
  message("raylib: USE_ROS enabled")
  add_definitions(-DUSE_ROS)
ENDIF(USE_ROS)

file(GLOB_RECURSE SOURCE_FILES
	${CMAKE_SOURCE_DIR}/raylib/src/*.cpp
)

file(GLOB_RECURSE INCLUDE_FILES
	${CMAKE_SOURCE_DIR}/raylib/include/*.h
)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  message_generation
  std_msgs
)
find_package(Eigen3 REQUIRED)
if (USE_LAS)
  find_package(libLAS REQUIRED)
endif()
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})

add_library(raylib ${SOURCE_FILES})
set_property(TARGET raylib PROPERTY POSITION_INDEPENDENT_CODE ON)

if (BUILD_SHARED_LIBS)
	target_compile_definitions(raylib PRIVATE "COMPILE_AS_DLL")
endif()

target_include_directories(raylib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/raylib/include>
  $<INSTALL_INTERFACE:include/raylib>
  ${EIGEN3_INCLUDE_DIR}
  ${catkin_INCLUDE_DIRS}
)
set_target_properties(raylib PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib
)
target_link_libraries(raylib ${catkin_LIBRARIES})
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
