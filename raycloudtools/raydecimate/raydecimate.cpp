// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Decimate a ray cloud spatially or temporally" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << std::endl;
  std::cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::IntArgument num_rays(1,100);
  ray::DoubleArgument vox_width(0.01, 100.0);
  ray::ValueKeyChoice quantity({&vox_width, &num_rays}, {"cm", "rays"});
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &quantity}))
    usage();

  std::string file = argv[1];
  ray::Cloud cloud;
  if (!cloud.load(file))
    usage();

  ray::Cloud new_cloud;
  std::string type = argv[3];
  if (type == "cm")
  {
    cloud.decimate(0.01 * std::stod(argv[2]));
  }
  else if (type == "rays")
  {
    int decimation = std::stoi(argv[2]);
    for (int i = 0; i < (int)cloud.ends.size(); i += decimation)
      new_cloud.addRay(cloud, i);
    cloud = new_cloud;
  }
  else
    usage(false);
  std::string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);
  cloud.save(file_stub + "_decimated.ply");
  return true;
}
