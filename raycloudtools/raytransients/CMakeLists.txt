set(SOURCES
  raytransients.cpp
)

ras_add_executable(raytransients
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
