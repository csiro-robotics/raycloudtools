// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raydebugdraw.h"

void usage(int exit_code = 1)
{
  std::cout << "Convert a point cloud and trajectory file into a ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayconvert pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << std::endl;
  std::cout << "                                            trajectoryfile is a text file in time,x,y,z,qw,qx,qy,qz format"
       << std::endl;
  std::cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file, trajectory_file;
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &trajectory_file}))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name(), trajectory_file.name()))
    usage();

  std::string save_file = cloud_file.nameStub();
  if (cloud_file.nameExt() == "ply")
    save_file += "_raycloud";
  cloud.save(save_file + ".ply");

  return 0;
}
