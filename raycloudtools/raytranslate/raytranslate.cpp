// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raymesh.h"

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
  std::cout << "raytranslate raycloud subtract ground_mesh.ply  - translate vertically to remove ground_mesh heights" << std::endl;
  std::cout << "raytranslate raycloud add ground_mesh.ply- translate vertically to add ground_mesh heights" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayTranslate(int argc, char *argv[])
{
  ray::FileArgument cloud_file, ground_file;
  ray::TextArgument subtract("subtract"), add("add");
  ray::Vector3dArgument translation3;
  ray::Vector4dArgument translation4;

  bool vec3_format = ray::parseCommandLine(argc, argv, { &cloud_file, &translation3 });
  bool vec4_format = ray::parseCommandLine(argc, argv, { &cloud_file, &translation4 });
  bool ground_subtract_format = ray::parseCommandLine(argc, argv, { &cloud_file, &subtract, &ground_file });
  bool ground_add_format = ray::parseCommandLine(argc, argv, { &cloud_file, &add, &ground_file });
  if (!vec3_format && !vec4_format && !ground_subtract_format && !ground_add_format)
    usage();

  Eigen::Vector3d translation(0, 0, 0);
  double time_delta = 0.0;
  ray::Mesh ground_mesh;
  ray::Cloud::Info info;
  Eigen::ArrayXXd ground_heights;  
  double voxel_width = 0.25;
  if (vec3_format)
    translation = translation3.value();
  else if (vec4_format)
  {
    translation = translation4.value().head<3>();
    time_delta = translation4.value()[3];
  }
  else
  {
    if (!ray::readPlyMesh(ground_file.name(), ground_mesh))
    {
      usage();
    }
    if (!ray::Cloud::getInfo(cloud_file.name(), info))
    {
      usage();
    }
    ground_mesh.toHeightField(ground_heights, info.ends_bound.min_bound_, info.ends_bound.max_bound_, voxel_width);
  }

  const std::string temp_name = cloud_file.nameStub() + "~.ply";  // tilde is a common suffix for temporary files

  auto translate = [&](Eigen::Vector3d &start, Eigen::Vector3d &end, double &time, ray::RGBA &) 
  {
    if (ground_subtract_format || ground_add_format)
    {
      Eigen::Vector3i index = ((end - info.ends_bound.min_bound_) / voxel_width).cast<int>();
      double height = ground_heights(index[0], index[1]);
      if (ground_subtract_format)
      {
        start[2] -= height;
        end[2] -= height;
      }
      else 
      {
        start[2] += height;
        end[2] += height;        
      }
    }
    else
    {
      start += translation;
      end += translation;
      time += time_delta;
    }
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