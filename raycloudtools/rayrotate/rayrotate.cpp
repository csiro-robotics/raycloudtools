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
  std::cout << "Rotate a raycloud about the origin" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrotate raycloud 30,0,0  - rotation (rx,ry,rz) is a rotation vector in degrees:" << std::endl;
  std::cout << "                             so this example rotates the cloud by 30 degrees in the x axis." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayRotate(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument rotation_arg(-360, 360);
  if (!ray::parseCommandLine(argc, argv, { &cloud_file, &rotation_arg }))
    usage();

  Eigen::Vector3d rot = rotation_arg.value();
  const double angle = rot.norm();
  rot /= angle;
  Eigen::Quaterniond rotation(Eigen::AngleAxisd(angle * ray::kPi / 180.0, rot));

  const std::string temp_name = cloud_file.name() + "~";  // tilde is a common suffix for temporary files

  auto rotate = [&](Eigen::Vector3d &start, Eigen::Vector3d &end, double &, ray::RGBA &) {
    start = rotation * start;
    end = rotation * end;
  };
  if (!ray::convertCloud(cloud_file.name(), temp_name, rotate))
    usage();

  std::rename(temp_name.c_str(), cloud_file.name().c_str());
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayRotate, argc, argv);
}