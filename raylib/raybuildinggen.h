// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBUILDINGGEN_H
#define RAYLIB_RAYBUILDINGGEN_H

#include "raylib/raycuboid.h"
#include "raylib/raylibconfig.h"
#include "rayutils.h"

namespace ray
{
/// Building raycloud generation class. Generates the attributes of a ray cloud for a randomly generated
/// building, containing doors, windows, boxes and cupboards
class RAYLIB_EXPORT BuildingGen
{
public:
  /// randomly generated building. The random seed can be chosen using @c srand()
  void generate();

  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }
  inline const std::vector<bool> rayBounded() const { return ray_bounded_; }

  struct RAYLIB_EXPORT BuildingParams
  {
    inline BuildingParams()
    {
      room_scales = Eigen::Vector3d(5.0, 5.0, 2.7);
      door_height = 2.0;
      door_width = 0.5;
      wall_width = 0.2;
      floor_width = 0.5;
      distinct_floor_likelihood = 0.55;
      outer_wall_width = 0.3;
      window_width = 1.2;
      window_height = 1.5;
      table_density = 0;     // initialised later
      cupboard_density = 0;  // initialised later
    }
    Eigen::Vector3d room_scales;
    double door_height;
    double door_width;
    double wall_width;
    double floor_width;
    double distinct_floor_likelihood;  // high has separated floors, low is a mixture of floor heights
    double outer_wall_width;
    double window_width;
    double window_height;
    double table_density;     // items per square metre
    double cupboard_density;  // items per metre along wall
  };

  inline BuildingParams &buildingParameters() { return params_; }
  inline const BuildingParams &buildingParameters() const { return params_; }

private:
  void splitRoom(const Cuboid &cuboid, std::vector<Cuboid> &cuboids, std::vector<std::vector<Cuboid>> &furniture);
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
  std::vector<bool> ray_bounded_;
  BuildingParams params_;
};
}  // namespace ray

#endif  // RAYLIB_RAYBUILDINGGEN_H
