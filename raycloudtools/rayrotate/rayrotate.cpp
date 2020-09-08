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
  std::cout << "Rotate a raycloud about the origin" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrotate raycloud 30,0,0  - rotation (rx,ry,rz) is a rotation vector in degrees:" << std::endl;
  std::cout << "                             so this example rotates the cloud by 30 degrees in the x axis." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument rotation;
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &rotation}))
    usage();
    
  ray::Pose pose;
  pose.position = Eigen::Vector3d(0, 0, 0);

  double angle = rotation.value.norm();
  rotation.value /= angle;
  pose.rotation = Eigen::Quaterniond(Eigen::AngleAxisd(angle * ray::kPi / 180.0, rotation.value));

  ray::Cloud cloud;
  cloud.load(cloud_file.name);
  cloud.transform(pose, 0.0);
  cloud.save(cloud_file.name);
  return true;
}
