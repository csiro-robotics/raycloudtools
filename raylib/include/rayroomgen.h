#pragma once
#include "rayutils.h"

namespace Ray
{
struct RoomGen
{
  void generate();
  vector<Eigen::Vector3d> rayStarts, rayEnds;
};
}
