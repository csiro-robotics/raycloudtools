set(SOURCES
  rayextract.cpp
)
add_compile_options(-ggdb -O0)

ras_add_executable(rayextract
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
