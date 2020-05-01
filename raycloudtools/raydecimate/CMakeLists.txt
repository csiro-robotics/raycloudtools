set(SOURCES
  raydecimate.cpp
)

ras_add_executable(raydecimate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
