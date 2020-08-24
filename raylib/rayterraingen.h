// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTERRAINGEN_H
#define RAYLIB_RAYTERRAINGEN_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

namespace ray
{
/// Terrain raycloud generation class. Generates the ray cloud attributes for a a random section of hilly terrain, 
/// as though scanned from a circular path over the terrain
class RAYLIB_EXPORT TerrainGen
{
public:
  /// terrain generation function. The random seed can be specified with @c srand()
  void generate();
  void generateRoughGround(double width, double roughness);

  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }
private:
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
};
}  // namespace ray

#endif  // RAYLIB_RAYTERRAINGEN_H
