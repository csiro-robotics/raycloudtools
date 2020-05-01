set(SOURCES
  raycolour.cpp
)

ras_add_executable(raycolour
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
