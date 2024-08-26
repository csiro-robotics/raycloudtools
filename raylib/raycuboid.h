// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCUBOID_H
#define RAYLIB_RAYCUBOID_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

namespace ray
{
/// Class for intersection tests on axis-aligned cuboids. These are closed intervals, their boundary is
/// included on intersection.
class RAYLIB_EXPORT Cuboid
{
public:
  Cuboid() {}
  Cuboid(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound);

  /// if the ray defined by @c start, @c dir and @c depth intersects the cuboid, then
  /// @c depth is updated with the intersection depth, and the function returns true
  /// @c positive_box true has outward normals and intersects with the near box faces
  /// @c positive_box false has inward normals and intersects with the far box faces
  bool intersectsRay(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth, bool positive_box) const;

  /// a negative box has inwards normals, so the ray intersects the back wall, filling in the @c depth argument
  bool intersects(const Eigen::Vector3d &pos) const;
  /// overlap of this cuboid with another, 'kissing' cuboids count as intersection.
  inline bool overlaps(const Cuboid &other) const
  {
    bool outside = other.min_bound_[0] > max_bound_[0] || other.min_bound_[1] > max_bound_[1] ||
                   other.min_bound_[2] > max_bound_[2] || other.max_bound_[0] < min_bound_[0] ||
                   other.max_bound_[1] < min_bound_[1] || other.max_bound_[2] < min_bound_[2];
    return !outside;
  }

  /// clip ray to cuboid. Return false if no ray left.
  bool clipRay(Eigen::Vector3d &start, Eigen::Vector3d &end, double eps = 0.0) const;

  Eigen::Vector3d min_bound_, max_bound_;
};
}  // namespace ray

#endif  // RAYLIB_RAYCUBOID_H
