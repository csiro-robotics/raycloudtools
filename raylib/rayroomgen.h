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
class RAYLIB_EXPORT RoomGen
{
public:
  void generate();
  std::vector<Eigen::Vector3d> ray_starts, ray_ends;
  std::vector<bool> ray_bounded;
};
}  // namespace ray

#endif  // RAYLIB_RAYROOMGEN_H
