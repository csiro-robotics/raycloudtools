// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raymesh.h"
#include <set>
#if defined(USE_QHULL)

namespace RAY
{
struct ConvexHull
{
  ConvexHull(const std::vector<Eigen::Vector3d> &points);

  void growOutwards(double maxCurvature);
  void growInwards(double maxCurvature);
  void growInDirection(double maxCurvature, const Eigen::Vector3d &dir);
  void growUpwards(double maxCurvature){ growInDirection(maxCurvature, Eigen::Vector3d(0,0,1)); }
  void growTopDown(double maxCurvature){ growInDirection(maxCurvature, Eigen::Vector3d(0,0,-1)); }

  Mesh mesh;
 
protected:
  void construct(const std::vector<Eigen::Vector3d> &points, const Eigen::Vector3d ignoreDirection);
};
}
#endif