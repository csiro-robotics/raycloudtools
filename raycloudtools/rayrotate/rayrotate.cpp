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
  std::cout << "Rotate a raycloud about the origin" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrotate raycloud 30,0,0  - rotation (rx,ry,rz) is a rotation vector in degrees:" << std::endl;
  std::cout << "                             so this example rotates the cloud by 30 degrees in the x axis." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();

  std::string file = argv[1];
  ray::Pose pose;
  pose.position = Eigen::Vector3d(0, 0, 0);

  std::stringstream ss(argv[2]);
  Eigen::Vector3d axis;
  ss >> axis[0];
  ss.ignore(1);
  ss >> axis[1];
  ss.ignore(1);
  ss >> axis[2];

  double angle = axis.norm();
  axis /= angle;
  pose.rotation = Eigen::Quaterniond(Eigen::AngleAxisd(angle * ray::kPi / 180.0, axis));

  ray::Cloud cloud;
  cloud.load(file);
  cloud.transform(pose, 0.0);
  cloud.save(file);
  return true;
}
