// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYROOMGEN_H
#define RAYLIB_RAYROOMGEN_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

namespace ray
{
struct RAYLIB_EXPORT RoomGen
{
  void generate();
  std::vector<Eigen::Vector3d> ray_starts, ray_ends;
  std::vector<bool> ray_bounded;
};
}

#endif // RAYLIB_RAYROOMGEN_H
