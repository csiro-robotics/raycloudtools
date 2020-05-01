set(SOURCES
  raycreate.cpp
)

ras_add_executable(raycreate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
