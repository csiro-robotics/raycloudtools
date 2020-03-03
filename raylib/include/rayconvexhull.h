#pragma once
#include "rayutils.h"
#include <set>

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

  std::vector<Eigen::Vector3d> points;
  std::vector<Eigen::Vector3d> vertices;
protected:
  void construct(const std::vector<Eigen::Vector3d> &points, const Eigen::Vector3d ignoreDirection);
};
}