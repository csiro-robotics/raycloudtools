set(SOURCES
  raywrap.cpp
)

ras_add_executable(raywrap
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
