set(SOURCES
  rayexport.cpp
)

ras_add_executable(rayexport
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
