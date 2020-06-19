set(SOURCES
  rayrestore.cpp
)

ras_add_executable(rayrestore
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
