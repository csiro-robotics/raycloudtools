set(SOURCES
  raysmooth.cpp
)

ras_add_executable(raysmooth
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
