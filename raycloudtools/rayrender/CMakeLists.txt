set(SOURCES
  rayrender.cpp
)
ras_add_executable(rayrender
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
