// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYMESH_H
#define RAYLIB_RAYMESH_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"

namespace ray
{

struct RAYLIB_EXPORT Mesh
{
  std::vector<Eigen::Vector3d> vertices;
  std::vector<Eigen::Vector3i> index_list;

  void splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside);
};


}

#endif // RAYLIB_RAYMESH_H
