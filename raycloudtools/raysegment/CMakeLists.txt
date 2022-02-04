set(SOURCES
  raysegment.cpp
)

ras_add_executable(raysegment
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
