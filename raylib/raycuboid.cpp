// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycuboid.h"
namespace ray
{

Cuboid::Cuboid(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound)
{
  min_bound_ = min_bound;
  max_bound_ = max_bound;
}

bool Cuboid::rayIntersectBox(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth) const
{
  double max_near_d = 0;
  double min_far_d = std::numeric_limits<double>::max();
  Eigen::Vector3d centre = (min_bound_ + max_bound_) / 2.0;
  Eigen::Vector3d extent = (max_bound_ - min_bound_) / 2.0;
  Eigen::Vector3d to_centre = centre - start;
  for (int ax = 0; ax < 3; ax++)
  {
    double s = dir[ax] > 0.0 ? 1.0 : -1.0;
    double near_d = (to_centre[ax] - s * extent[ax]) / dir[ax];
    double far_d = (to_centre[ax] + s * extent[ax]) / dir[ax];

    max_near_d = std::max(max_near_d, near_d);
    min_far_d = std::min(min_far_d, far_d);
  }
  if (max_near_d > 0.0 && max_near_d < depth && max_near_d < min_far_d)
  {
    depth = max_near_d;
    return true;
  }
  return false;
}

bool Cuboid::rayIntersectNegativeBox(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth) const
{
  double max_near_d = 0;
  double min_far_d = std::numeric_limits<double>::max();
  Eigen::Vector3d centre = (min_bound_ + max_bound_) / 2.0;
  Eigen::Vector3d extent = (max_bound_ - min_bound_) / 2.0;
  Eigen::Vector3d to_centre = centre - start;
  for (int ax = 0; ax < 3; ax++)
  {
    double s = dir[ax] > 0.0 ? 1.0 : -1.0;
    double near_d = (to_centre[ax] - s * extent[ax]) / dir[ax];
    double far_d = (to_centre[ax] + s * extent[ax]) / dir[ax];

    max_near_d = std::max(max_near_d, near_d);
    min_far_d = std::min(min_far_d, far_d);
  }
  if (max_near_d < min_far_d && min_far_d < depth)
  {
    depth = min_far_d;
    return true;
  }
  return false;
}

bool Cuboid::intersects(const Eigen::Vector3d &pos) const
{
  return pos[0] >= min_bound_[0] && pos[1] >= min_bound_[1] && pos[2] >= min_bound_[2] && pos[0] <= max_bound_[0] &&
          pos[1] <= max_bound_[1] && pos[2] <= max_bound_[2];
}

bool Cuboid::overlaps(const Cuboid &other) const
{
  bool outside = other.min_bound_[0] > max_bound_[0] || other.min_bound_[1] > max_bound_[1] || other.min_bound_[2] > max_bound_[2] ||
                  other.max_bound_[0] < min_bound_[0] || other.max_bound_[1] < min_bound_[1] || other.max_bound_[2] < min_bound_[2];
  return !outside;
}

} // ray