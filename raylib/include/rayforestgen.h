// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raytreegen.h"

namespace RAY
{
struct ForestGen
{
  void make(double randomFactor = 0.0);
  void generateRays(double rayDensity);
  std::vector<Eigen::Vector3d> getCanopy();
  std::vector<Eigen::Vector3d> getPointCloud();

  std::vector<TreeGen> trees;
};
}
