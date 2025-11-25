set(SOURCES
  raydiff.cpp
)
ras_add_executable(raydiff
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
