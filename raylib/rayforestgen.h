// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFORESTGEN_H
#define RAYLIB_RAYFORESTGEN_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raytreegen.h"

namespace ray
{
struct RAYLIB_EXPORT ForestGen
{
  void make(double random_factorr = 0.0);
  void generateRays(double ray_densityy);
  std::vector<Eigen::Vector3d> getCanopy();
  std::vector<Eigen::Vector3d> getPointCloud();

  std::vector<TreeGen> trees;
};
}  // namespace ray

#endif  // RAYLIB_RAYFORESTGEN_H
