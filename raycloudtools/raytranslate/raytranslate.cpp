// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Translate a raycloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << std::endl;
  std::cout << "                      0,0,1,24.3 - optional 4th component translates time" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayTranslate(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument translation3;
  ray::Vector4dArgument translation4;

  bool vec3 = ray::parseCommandLine(argc, argv, { &cloud_file, &translation3 });
  bool vec4 = ray::parseCommandLine(argc, argv, { &cloud_file, &translation4 });
  if (!vec3 && !vec4)
    usage();

  Eigen::Vector3d translation(0, 0, 0);
  double time_delta = 0.0;
  if (vec3)
    translation = translation3.value();
  else
  {
    translation = translation4.value().head<3>();
    time_delta = translation4.value()[3];
  }

  const std::string temp_name = cloud_file.nameStub() + "~.ply";  // tilde is a common suffix for temporary files

  auto translate = [&](Eigen::Vector3d &start, Eigen::Vector3d &end, double &time, ray::RGBA &) {
    start += translation;
    end += translation;
    time += time_delta;
  };
  if (!ray::convertCloud(cloud_file.name(), temp_name, translate))
    usage();

  std::rename(temp_name.c_str(), cloud_file.name().c_str());

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayTranslate, argc, argv);
}