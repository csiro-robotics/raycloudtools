set(SOURCES
  raysplit.cpp
)
add_compile_options(-ggdb -O0)
ras_add_executable(raysplit
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
