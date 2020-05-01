set(SOURCES
  raycombine.cpp
)

ras_add_executable(raycombine
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
