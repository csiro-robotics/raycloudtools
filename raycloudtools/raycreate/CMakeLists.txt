set(SOURCES
  raycreate.cpp
)
add_compile_options(-ggdb -O0)
ras_add_executable(raycreate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
