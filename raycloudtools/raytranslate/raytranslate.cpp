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
#define DETAILED_GROUND_SUBTRACTION // ground lookup per-point. Otherwise it quantises to a height grid

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Translate a raycloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << std::endl;
  std::cout << "                      0,0,1,24.3 - optional 4th component translates time" << std::endl;
  std::cout << "                      subtract ground_mesh.ply  - translate vertically to remove ground_mesh heights" << std::endl;
  std::cout << "                      add ground_mesh.ply- translate vertically to add ground_mesh heights" << std::endl;
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
  Eigen::Vector3d min_bound;
  Eigen::Vector3d max_bound;
  #if defined DETAILED_GROUND_SUBTRACTION
  std::vector<ray::Triangle> triangles;
  ray::Grid<ray::Triangle *> grid;
  #else
  Eigen::ArrayXXd ground_heights;  
  #endif
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
    min_bound = info.ends_bound.min_bound_;
    max_bound = info.ends_bound.max_bound_;
    // the toHeightfield function only casts rays from max_bound[2] down to min_bound[2], so set these according to
    // ground_file's extents, not cloud_file's extents.
    min_bound[2] = std::numeric_limits<double>::max();
    max_bound[2] = std::numeric_limits<double>::lowest();
    for (auto &vert: ground_mesh.vertices())
    {
      min_bound[2] = std::min(min_bound[2], vert[2]);
      max_bound[2] = std::max(max_bound[2], vert[2]);
    }
    min_bound[2] -= 0.1;
    max_bound[2] += 0.1;
    #if defined DETAILED_GROUND_SUBTRACTION
    Eigen::Vector3d box_min = min_bound;
    Eigen::Vector3d box_max = max_bound;
    box_max[2] = box_min[2] + 0.5 * voxel_width;  // ensure that the grid is only 1 voxel high
    triangles.resize(ground_mesh.indexList().size());
    for (int i = 0; i < (int)ground_mesh.indexList().size(); i++)
    {
      ray::Triangle &tri = triangles[i];
      for (int j = 0; j < 3; j++) 
      {
        tri.corners[j] = ground_mesh.vertices()[ground_mesh.indexList()[i][j]];
      }
      tri.tested = false;
      tri.normal = (tri.corners[1] - tri.corners[0]).cross(tri.corners[2] - tri.corners[0]);
    }

    // put the triangles into a grid
    grid.init(box_min, box_max, voxel_width);
    for (auto &tri : triangles)
    {
      Eigen::Vector3d tri_min = (ray::minVector(tri.corners[0], ray::minVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
      Eigen::Vector3d tri_max = (ray::maxVector(tri.corners[0], ray::maxVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
      for (int x = (int)tri_min[0]; x <= (int)tri_max[0]; x++)
      {
        for (int y = (int)tri_min[1]; y <= (int)tri_max[1]; y++) 
        {
          grid.insert(x, y, 0, &tri);
        }
      }
    }
    #else
    ground_mesh.toHeightField(ground_heights, min_bound, max_bound, voxel_width);
    #endif
  }

  const std::string temp_name = cloud_file.nameStub() + "~.ply";  // tilde is a common suffix for temporary files
  int num_missed_triangles = 0, num_totally_missed = 0;

  auto translate = [&](Eigen::Vector3d &start, Eigen::Vector3d &end, double &time, ray::RGBA &) 
  {
    if (ground_subtract_format || ground_add_format)
    {
      #if defined DETAILED_GROUND_SUBTRACTION
      Eigen::Vector3i index = ((end - info.ends_bound.min_bound_) / voxel_width).cast<int>();
      int x = index[0];
      int y = index[1];

      // now look up the triangle for each pixel centre
      Eigen::Vector3d pos_top = min_bound + voxel_width * (Eigen::Vector3d((double)x + 0.5, (double)y + 0.5, 0));
      Eigen::Vector3d pos_base = pos_top;
      pos_top[2] = max_bound[2];
      pos_base[2] = min_bound[2];
      auto &tris = grid.cell(x, y, 0).data;
      // search the triangles in this cell 'bucket'
      double height = 0.0;
      double mean_depth = 0.0;
      double mean_count = 0.0;
      for (auto &tri : tris)
      {
        double depth = 0.0;
        if (tri->intersectsRay(pos_top, pos_base, depth))
        {
          height = pos_top[2] + (pos_base[2] - pos_top[2]) * depth;
          break;
        }
        if (depth != 0.0)
        {
          mean_depth += depth;
          mean_count++;
        }
      }
      if (height == 0.0)
      {
        if (mean_count > 0)
        {
          mean_depth /= mean_count;
          height = pos_top[2] + (pos_base[2] - pos_top[2]) * mean_depth;
          num_missed_triangles++;
        }
        else 
        {
          num_totally_missed++;
        }
      }
      #else
      Eigen::Vector3i index = ((end - info.ends_bound.min_bound_) / voxel_width).cast<int>();
      double height = ground_heights(index[0], index[1]);
      #endif
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

  if (num_missed_triangles > 2)
  {
    std::cout << num_missed_triangles << " cloud points not overlapping triangles, if this is large the ground file may not be a lateral coverage of the ray cloud" << std::endl;
  }
  if (num_totally_missed > 0)
  {
    std::cout << "Warning: " << num_totally_missed << " points have no laterally overlapping triangles, so the ground_file is not a full coverage. Point translation ignored" << std::endl;
  }
  std::rename(temp_name.c_str(), cloud_file.name().c_str());

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayTranslate, argc, argv);
}