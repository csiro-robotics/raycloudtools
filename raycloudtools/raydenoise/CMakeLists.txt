set(SOURCES
  raydenoise.cpp
)

ras_add_executable(raydenoise
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
