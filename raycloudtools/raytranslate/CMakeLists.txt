set(SOURCES
  raytranslate.cpp
)
add_compile_options(-ggdb -O0)
ras_add_executable(raytranslate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
