#pragma once
#include "rayutils.h"

namespace RAY
{
struct RoomGen
{
  void generate();
  std::vector<Eigen::Vector3d> rayStarts, rayEnds;
};
}
