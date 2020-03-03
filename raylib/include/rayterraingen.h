#pragma once
#include "rayutils.h"

namespace Ray
{
struct TerrainGen
{
  void generate();
  std::vector<Eigen::Vector3d> rayStarts, rayEnds;
};
}
