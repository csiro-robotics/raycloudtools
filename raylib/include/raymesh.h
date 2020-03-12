// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raycloud.h"

namespace RAY
{

struct Mesh
{
  std::vector<Eigen::Vector3d> vertices;
  std::vector<Eigen::Vector3i> indexList;

  void splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside);
};


}
