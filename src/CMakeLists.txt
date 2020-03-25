add_compile_options(-g -Wall -O3 -fopenmp)

OPTION(USE_LAS
  "Build the project using LAS"
  OFF)
IF(USE_LAS)
  message("LAS enabled")
  add_definitions(-DUSE_LAS)
ENDIF(USE_LAS)

OPTION(USE_QHULL
"Build the project using QHull"
OFF)
IF(USE_QHULL)
  message("qhull enabled")
  add_definitions(-DUSE_QHULL)
ENDIF(USE_QHULL)


add_definitions(-DNDEBUG )
find_package(catkin REQUIRED COMPONENTS
roscpp
message_generation
std_msgs
)
if (USE_LAS)
  find_package(libLAS REQUIRED)
endif()

find_package(libnabo REQUIRED)
include_directories(../raylib/include ${catkin_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})

# qhull
if (USE_QHULL)
  find_path(TRAVERSIBILITY_QHULL_INCLUDE_DIR libqhullcpp/Qhull.h)
  if (NOT TRAVERSIBILITY_QHULL_INCLUDE_DIR)
    message(FATAL_ERROR "Unable to find libqhullcpp/Qhull.h")
  endif()
  find_library(TRAVERSIBILITY_LIBQHULLCPP qhullcpp libqhullcpp)
  if (NOT TRAVERSIBILITY_LIBQHULLCPP)
    message(FATAL_ERROR "You need to have the qhull cpp library installed.")
  endif()
  SET(qhullDir ../qhull-2019.1)
  INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
  INCLUDE_DIRECTORIES(${qhullDir}/src)
  LINK_DIRECTORIES(${qhullDir}/build)
  INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
  INCLUDE_DIRECTORIES(include)
  SET(qhullLibs qhullcpp qhull_r)
endif()

set(TARGETS
  rayalign
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
  raycombine
)
foreach(TARGET ${TARGETS})
	add_executable(${TARGET} ${TARGET}.cpp)
endforeach(TARGET)

foreach(TARGET ${TARGETS})
	target_link_libraries(${TARGET} raylib ${catkin_LIBRARIES} ${libLAS_LIBRARIES} ${libnabo_LIBRARIES} ${qhullLibs})
	set_target_properties(${TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../bin)
	set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${TARGET})
	install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
endforeach(TARGET)
