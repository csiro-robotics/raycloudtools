set(SOURCES
  raysplit.cpp
)

ras_add_executable(raysplit
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
