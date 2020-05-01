set(SOURCES
  rayconvert.cpp
)

ras_add_executable(rayconvert
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
