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
  std::cout << "Translate a raycloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << std::endl;
  std::cout << "                      0,0,1,24.3 - optional 4th component translates time" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();

  std::string file = argv[1];
  ray::Pose pose;
  std::stringstream ss(argv[2]);
  Eigen::Vector3d vec;
  ss >> vec[0];
  ss.ignore(1);
  ss >> vec[1];
  ss.ignore(1);
  ss >> vec[2];
  pose.position = vec;

  double time_delta = 0.0;
  if (ss.peek() == ',')
  {
    ss.ignore(1);
    ss >> time_delta;
  }
  pose.rotation = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);

  ray::Cloud cloud;
  cloud.load(file);
#if 0  // test bending
  Eigen::Vector3d bend(pose.position[0], pose.position[1], 0.0);
  Eigen::Vector3d side = Eigen::Vector3d(bend[1], -bend[0], 0).normalized();
  for (auto &end: cloud.ends)
    end += bend*sqr(end.dot(side));
#else
  cloud.transform(pose, time_delta);
#endif
  cloud.save(file);
  return true;
}
