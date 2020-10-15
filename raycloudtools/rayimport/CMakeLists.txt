set(SOURCES
  rayimport.cpp
)

ras_add_executable(rayimport
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
