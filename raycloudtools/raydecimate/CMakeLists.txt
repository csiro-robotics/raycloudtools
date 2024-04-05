set(SOURCES
  raydecimate.cpp
)
add_compile_options(-ggdb -O0)
ras_add_executable(raydecimate
  LIBS raylib
  SOURCES ${SOURCES}
  PROJECT_FOLDER "raycloudtools"
)
