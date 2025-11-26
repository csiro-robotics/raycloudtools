set(SOURCES
  raydiff.cpp
)
add_compile_options(-ggdb -O0)
ras_add_executable(raydiff
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
