// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYROOMGEN_H
#define RAYLIB_RAYROOMGEN_H

#include "raylib/raycuboid.h"
#include "raylib/raylibconfig.h"
#include "rayutils.h"

namespace ray
{
/// Room raycloud generation class. Generates the attributes of a ray cloud for a randomly generated
/// single room, containing a table, a wardrobe, a door and a window
class RAYLIB_EXPORT RoomGen
{
public:
  /// randomly generated room. The random seed can be chosen using @c srand()
  void generate();

  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }
  inline const std::vector<bool> rayBounded() const { return ray_bounded_; }

private:
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
  std::vector<bool> ray_bounded_;
};
}  // namespace ray

#endif  // RAYLIB_RAYROOMGEN_H
