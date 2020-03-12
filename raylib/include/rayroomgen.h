// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
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
