// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include <raylib/raycloud.h>

#include <cstdlib>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Load a ray cloud using raycloudtools::raylib" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "useraylib raycloud" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    usage();
  }

  ray::Cloud raycloud;
  raycloud.load(argv[1]);

  std::cout << "Cloud " << argv[1] << " loaded with " << raycloud.starts.size() << " positions and "
            << raycloud.ends.size() << "points" << std::endl;

  return 0;
}
