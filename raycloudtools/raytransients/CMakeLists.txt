set(SOURCES
  raytransients.cpp
)

ras_add_executable(raytrainsients
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
