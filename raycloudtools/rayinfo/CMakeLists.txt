set(SOURCES
  rayinfo.cpp
)

ras_add_executable(rayinfo
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
