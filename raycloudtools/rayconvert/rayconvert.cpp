// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Convert a point cloud and trajectory file into a ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayconvert pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << std::endl;
  std::cout << "                                            trajectoryfile is a text file in time,x,y,z,qw,qx,qy,qz format"
       << std::endl;
  std::cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << std::endl;
  exit(exit_code);
}

#include "raylib/raydebugdraw.h"

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();

  std::string point_cloud_file = argv[1];
  std::string trajectory_file = argv[2];
  ray::Cloud cloud;
  cloud.load(point_cloud_file, trajectory_file);

  std::string save_file = point_cloud_file.substr(0, point_cloud_file.size() - 4);
  if (point_cloud_file.substr(point_cloud_file.size() - 4) == ".ply")
    save_file += "_raycloud";
  cloud.save(save_file + ".ply");

  return true;
}
