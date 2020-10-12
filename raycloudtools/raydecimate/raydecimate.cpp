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

void usage(int exit_code = 1)
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

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();

  ray::Cloud new_cloud;
  std::string type = quantity.selectedKey();
  if (type == "cm")
  {
    cloud.decimate(0.01 * vox_width.value());
  }
  else if (type == "rays")
  {
    int decimation = num_rays.value();
    for (int i = 0; i < (int)cloud.ends.size(); i += decimation)
      new_cloud.addRay(cloud, i);
    cloud = new_cloud;
  }
  else
    usage(false);
  cloud.save(cloud_file.nameStub() + "_decimated.ply");
  return 0;
}
