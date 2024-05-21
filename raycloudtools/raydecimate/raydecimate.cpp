// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raydecimation.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Decimate a ray cloud spatially or temporally" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm. A spatially even subsampling" << std::endl;
  std::cout << "raydecimate raycloud 4 rays - reduces to every fourth ray. A temporally even subsampling (if rays are chronological)" << std::endl;
  std::cout << "advanced methods not supported in rayrestore:" << std::endl;
  std::cout << "raydecimate raycloud 20 cm 64 rays - A maximum of 64 rays per cubic 20 cm. Retains small-scale details compared to spatial decimation" << std::endl;
  std::cout << "raydecimate raycloud 3 cm/m - reduces to ray ends spaced 3 cm apart for each metre of their length" << std::endl;
  // clang-format off
  exit(exit_code);
}

class Vector4iLess
{
public:
  bool operator()(const Eigen::Vector4i &a, const Eigen::Vector4i &b) const
  {
    if (a[3] != b[3])
      return a[3] < b[3];
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    return a[2] < b[2];
  }
};

// Decimates the ray cloud, spatially or in time
int rayDecimate(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::IntArgument num_rays(1, 10000);
  ray::DoubleArgument vox_width(0.01, 100.0);
  ray::DoubleArgument radius_per_length(0.01, 100.0);
  ray::ValueKeyChoice quantity({ &vox_width, &num_rays, &radius_per_length }, { "cm", "rays", "cm/m" });
  ray::TextArgument cm("cm"), rays("rays"); 
  bool standard_format = ray::parseCommandLine(argc, argv, { &cloud_file, &quantity });
  bool double_format = ray::parseCommandLine(argc, argv, { &cloud_file, &vox_width, &cm, &num_rays, &rays });
  if (!standard_format && !double_format)
    usage();

  bool res = false;
  if (double_format)
  {
    res = ray::decimateSpatioTemporal(cloud_file.nameStub(), vox_width.value(), num_rays.value());
  }
  else if (quantity.selectedKey() == "cm")
  {
    res = ray::decimateSpatial(cloud_file.nameStub(), vox_width.value());
  }
  else if (quantity.selectedKey() == "rays")
  {
    res = ray::decimateTemporal(cloud_file.nameStub(), num_rays.value());
  }
  else if (quantity.selectedKey() == "cm/m")
  {
    res = ray::decimateAngular(cloud_file.nameStub(), radius_per_length.value());
  }
  if (!res)
    usage();

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDecimate, argc, argv);
}