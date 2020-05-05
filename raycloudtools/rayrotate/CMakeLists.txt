set(SOURCES
  rayrotate.cpp
)

ras_add_executable(rayrotate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
