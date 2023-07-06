// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tom Lowe, Kazys Stepanas
#ifndef RAYELLIPSOID_H
#define RAYELLIPSOID_H

#include "raylib/raylibconfig.h"

#include <Eigen/Dense>

#include <cmath>
#include <vector>

namespace ray
{
class Cloud;
class Progress;

enum class RAYLIB_EXPORT IntersectResult
{
  Miss,
  Passthrough,
  Hit,
};

class RAYLIB_EXPORT Ellipsoid
{
public:
  Eigen::Vector3d pos;
  Eigen::Matrix3f eigen_mat;  // each row is a scaled eigenvector
  Eigen::Vector3f extents;
  double time;
  float opacity;  ///< A representation of certainty of this ellipsoid.
  size_t num_rays;
  size_t num_gone;
  bool transient;

  void clear();
  void setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals);

  IntersectResult intersect(const Eigen::Vector3d &start, const Eigen::Vector3d &end) const;
};

/// Convert the cloud into a list of ellipsoids, which represent a volume around each cloud point,
/// shaped by the distribution of its neighbouring points.
void RAYLIB_EXPORT generateEllipsoids(std::vector<Ellipsoid> *ellipsoids, Eigen::Vector3d *bounds_min,
                                      Eigen::Vector3d *bounds_max, const Cloud &cloud, Progress *progress = nullptr);

inline void Ellipsoid::clear()
{
  pos = Eigen::Vector3d::Zero();
  eigen_mat = Eigen::Matrix3f::Identity();
  extents = Eigen::Vector3f::Zero();
  time = 0.0;
  opacity = 0.0f;
  num_rays = num_gone = 0;
  transient = false;
}

inline void Ellipsoid::setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals)
{
  // This is approximate (slightly larger than minimal bounds), but
  // an exact bounding box is most likely non-analytic, and expensive to compute
  double max_rr = std::max(vals[0], std::max(vals[1], vals[2]));
  const Eigen::Vector3d &x = vecs.col(0);
  const Eigen::Vector3d &y = vecs.col(1);
  const Eigen::Vector3d &z = vecs.col(2);
  extents[0] = static_cast<float>(std::min(max_rr, std::abs(x[0]) * vals[0] + std::abs(y[0]) * vals[1] + std::abs(z[0]) * vals[2]));
  extents[1] = static_cast<float>(std::min(max_rr, std::abs(x[1]) * vals[0] + std::abs(y[1]) * vals[1] + std::abs(z[1]) * vals[2]));
  extents[2] = static_cast<float>(std::min(max_rr, std::abs(x[2]) * vals[0] + std::abs(y[2]) * vals[1] + std::abs(z[2]) * vals[2]));
}

inline IntersectResult Ellipsoid::intersect(const Eigen::Vector3d &start, const Eigen::Vector3d &end) const
{
  Eigen::Matrix3d eigen_matd = eigen_mat.cast<double>();
  const Eigen::Vector3d dir = end - start;
  // ray-ellipsoid intersection
  const Eigen::Vector3d to_sphere = pos - start;
  const Eigen::Vector3d ray = eigen_matd * dir;
  const double ray_length_sqr = ray.squaredNorm();
  const Eigen::Vector3d to = eigen_matd * to_sphere;

  double d = to.dot(ray) / ray_length_sqr;
  const double dist2 = (to - ray * d).squaredNorm();

  if (dist2 > 1.0)  // misses the ellipsoid
  {
    return IntersectResult::Miss;
  }

  const double along_dist = std::sqrt(1.0 - dist2);
  const double ray_length = std::sqrt(ray_length_sqr);
  d *= ray_length;
  if (ray_length < d - along_dist)  // doesn't reach the ellipsoid
  {
    return IntersectResult::Miss;
  }

  const double pass_distance = 0.05;
  double ratio = pass_distance / dir.norm();
  // last number requires rays to pass some way past the object
  const bool pass_through = ray_length * (1.0 - ratio) > d + along_dist;
  if (pass_through)
  {
    return IntersectResult::Passthrough;
  }
  return IntersectResult::Hit;
}

}  // namespace ray

#endif  // RAYELLIPSOID_H
