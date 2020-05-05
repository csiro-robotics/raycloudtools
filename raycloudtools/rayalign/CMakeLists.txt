set(SOURCES
  rayalign.cpp
)

ras_add_executable(rayalign
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
