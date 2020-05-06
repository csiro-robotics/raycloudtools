// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYINTERSECT_H
#define RAYINTERSECT_H

#include "raylib/raylibconfig.h"

#include <Eigen/Dense>

#include <algorithm>
#include <limits>

namespace ray
{
namespace intersect
{
/// Test if two axis aligned bounding boxes overlap.
inline bool aabbOverlap(const Eigen::Vector3d &aabb_min_a, const Eigen::Vector3d &aabb_max_a,
                        const Eigen::Vector3d &aabb_min_b, const Eigen::Vector3d &aabb_max_b, const double epsilon = 0)
{
  using Vector3b = Eigen::Matrix<bool, 3, 1>;
  const Vector3b max_less_min(aabb_max_a.x() + epsilon < aabb_min_b.x(),  //
                              aabb_max_a.y() + epsilon < aabb_min_b.y(),  //
                              aabb_max_a.z() + epsilon < aabb_min_b.z());
  const Vector3b min_greater_max(aabb_min_a.x() - epsilon < aabb_max_b.x(),  //
                                 aabb_min_a.y() - epsilon < aabb_max_b.y(),  //
                                 aabb_min_a.z() - epsilon < aabb_max_b.z());

  return !(max_less_min.x() || max_less_min.y() || max_less_min.z()) &&
         !(min_greater_max.x() || min_greater_max.y() || min_greater_max.z());
}

/// Test if @p point lies within the given axis aligned bounding box bounds.
inline bool pointInAabb(const Eigen::Vector3d &point, const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max,
                        const double epsilon = 0)
{
  return aabbOverlap(aabb_min, aabb_max, point, point, epsilon);
}

/// Helper function for @c rayAabb(). Calculates the hit time for a single axis.
inline double calcTimeVal(double limit, double origin, double direction);
/// Helper function for @c rayAabb(). Calculates the overlap between ranges @p a and @p b into @p overlap.
inline bool calcIntervalOverlap(const double a[2], const double b[2], double overlap[2]);
/// Helper function for @c rayAabb. Returns 1 if @p val is zero or positive, -1 if @p val is negative.
inline double sign(double val) { return val >= 0 ? 1.0f : -1.0f; }

/// Test a ray against an axis aligned bounding box, filling @p hit_times with the entry and exit times.
/// These times are ranges when @p ray is unit length, otheriwse they are scale factors. The return values
/// may be less than zero and greater than the length of @p ray indicating an intersection before the @p origin point
/// or beyond the @p ray length respectively.
inline bool rayAabb(const Eigen::Vector3d &origin, const Eigen::Vector3d &ray, const Eigen::Vector3d &aabb_min,
                    const Eigen::Vector3d &aabb_max, double hit_times[2])
{
  // Based on: https://tavianator.com/fast-branchless-raybounding-box-intersections/
  // Convert to ray format.
  const Eigen::Vector3d inv_dir(1.0 / ray.x(), 1.0 / ray.y(), 1.0 / ray.z());
  const Eigen::Vector3i sign(!!(inv_dir.x() < 0), !!(inv_dir.y() < 0), !!(inv_dir.z() < 0));
  const Eigen::Vector3d aabb[2] = { aabb_min, aabb_max };
  double tx[2], ty[2], tz[2];

  tx[0] = calcTimeVal(aabb[sign[0]].x(), origin.x(), inv_dir.x());
  tx[1] = calcTimeVal(aabb[1 - sign[0]].x(), origin.x(), inv_dir.x());
  ty[0] = calcTimeVal(aabb[sign[1]].y(), origin.y(), inv_dir.y());
  ty[1] = calcTimeVal(aabb[1 - sign[1]].y(), origin.y(), inv_dir.y());

  if (!calcIntervalOverlap(tx, ty, hit_times))
  {
    return false;
  }

  tz[0] = calcTimeVal(aabb[sign[2]].z(), origin.z(), inv_dir.z());
  tz[1] = calcTimeVal(aabb[1 - sign[2]].z(), origin.z(), inv_dir.z());

  if (!calcIntervalOverlap(hit_times, tz, hit_times))
  {
    return false;
  }

  return true;
  // bool intersected = false;
  // if (time_best[0] > 0)
  // {
  //   intersected = true;
  // }

  // if (time_best[1] < max_time)
  // {
  //   intersected = true;
  // }

  // return intersected;
}


inline double calcTimeVal(double limit, double origin, double direction)
{
  // Always performing the line on the return value nearly works, but occasionally
  // returns NaN instead of infinite results.
  if (direction == std::numeric_limits<double>::infinity())
  {
    return sign(limit - origin) * std::numeric_limits<double>::infinity();
  }
  if (direction == -std::numeric_limits<double>::infinity())
  {
    return sign(limit - origin) * -std::numeric_limits<double>::infinity();
  }
  return (limit - origin) * direction;
}


inline bool calcIntervalOverlap(const double a[2], const double b[2], double overlap[2])
{
  if (a[0] > b[1])
  {
    return false;
  }
  if (b[0] > a[1])
  {
    return false;
  }

  overlap[0] = std::max(a[0], b[0]);
  overlap[1] = std::min(a[1], b[1]);
  return true;
}
}  // namespace intersect
}  // namespace ray

#endif  // RAYINTERSECT_H
