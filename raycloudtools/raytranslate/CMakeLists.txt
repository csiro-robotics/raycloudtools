set(SOURCES
  raytranslate.cpp
)

ras_add_executable(raytranslate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
