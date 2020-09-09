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
  std::cout << "Translate a raycloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << std::endl;
  std::cout << "                      0,0,1,24.3 - optional 4th component translates time" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument translation3;
  ray::Vector4dArgument translation4;

  bool vec3 = ray::parseCommandLine(argc, argv, {&cloud_file, &translation3});
  bool vec4 = ray::parseCommandLine(argc, argv, {&cloud_file, &translation4});
  if (!vec3 && !vec4)
    usage();

  ray::Pose pose;
  double time_delta = 0.0;
  if (vec3)
    pose.position = translation3.value;
  else
  {
    pose.position = translation4.value.head<3>();
    time_delta = translation4.value[3];
  }
  pose.rotation = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name))
    usage();
#if 0  // test bending
  Eigen::Vector3d bend(pose.position[0], pose.position[1], 0.0);
  Eigen::Vector3d side = Eigen::Vector3d(bend[1], -bend[0], 0).normalized();
  for (auto &end: cloud.ends)
    end += bend*sqr(end.dot(side));
#else
  cloud.transform(pose, time_delta);
#endif
  cloud.save(cloud_file.name);
  return true;
}
