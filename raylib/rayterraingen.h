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
// parameters for specifying the form of terrain
struct TerrainParams
{
  TerrainParams()
    : point_density(400)
    , ray_height(2.0)
    , range_noise(0.03)
    , walk_radius(4.0)
  {}
  double point_density;  // in points per square metre
  double ray_height;
  double range_noise;
  double walk_radius;  // radius of path of trajectory
};

/// Terrain raycloud generation class. Generates the ray cloud attributes for a a random section of hilly terrain,
/// as though scanned from a circular path over the terrain
class RAYLIB_EXPORT TerrainGen
{
public:
  /// terrain generation function. The random seed can be specified with @c srand()
  void generate(const TerrainParams &params = TerrainParams());

  /// generate terrain from a mesh file @c filename, and according to the @c params
  bool generateFromFile(const std::string &filename, const TerrainParams &params = TerrainParams());

  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }

private:
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
};
}  // namespace ray

#endif  // RAYLIB_RAYTERRAINGEN_H
