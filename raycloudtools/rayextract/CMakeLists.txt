set(SOURCES
  rayextract.cpp
)

ras_add_executable(rayextract
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
