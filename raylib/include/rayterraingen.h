#pragma once
#include "rayutils.h"

namespace RAY
{
struct TerrainGen
{
  void generate();
  std::vector<Eigen::Vector3d> rayStarts, rayEnds;
};
}
