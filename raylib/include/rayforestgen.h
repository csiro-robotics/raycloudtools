#pragma once
#include "rayutils.h"
#include "raytreegen.h"

namespace Ray
{
struct ForestGen
{
  void make(double randomFactor = 0.0);
  void generateRays(double rayDensity);
  vector<Vector3d> getCanopy();
  vector<Vector3d> getPointCloud();

  vector<TreeGen> trees;
};
}
