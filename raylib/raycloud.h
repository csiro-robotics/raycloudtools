// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUD_H
#define RAYLIB_RAYCLOUD_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"
#include "raygrid.h"
#include <set>

namespace ray
{
class Progress;

typedef Eigen::Matrix<double, 6, 1> Vector6i;

struct RAYLIB_EXPORT Vector6iLess
{
  inline bool operator()(const Vector6i &a, const Vector6i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    if (a[2] != b[2])
      return a[2] < b[2];
    if (a[3] != b[3])
      return a[3] < b[3];
    if (a[4] != b[4])
      return a[4] < b[4];
    return a[5] < b[5];
  }
};

enum class RAYLIB_EXPORT IntersectResult
{
  Miss,
  Passthrough,
  Hit,
};

struct RAYLIB_EXPORT Ellipsoid
{
  Eigen::Vector3d pos;
  Eigen::Matrix3d eigen_mat;  // each row is a scaled eigenvector
  double time;
  Eigen::Vector3d extents;
  double opacity;
  double planarity;
  size_t num_rays;
  size_t num_gone;
  bool transient;

  inline void setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals)
  {
    // This is approximate (slightly larger than minimal bounds), but
    // an exact bounding box is most likely non-analytic, and expensive to compute
    double max_rr = std::max(vals[0], std::max(vals[1], vals[2]));
    const Eigen::Vector3d &x = vecs.col(0);
    const Eigen::Vector3d &y = vecs.col(1);
    const Eigen::Vector3d &z = vecs.col(2);
    extents[0] = std::min(max_rr, abs(x[0]) * vals[0] + abs(y[0]) * vals[1] + abs(z[0]) * vals[2]);
    extents[1] = std::min(max_rr, abs(x[1]) * vals[0] + abs(y[1]) * vals[1] + abs(z[1]) * vals[2]);
    extents[2] = std::min(max_rr, abs(x[2]) * vals[0] + abs(y[2]) * vals[1] + abs(z[2]) * vals[2]);
  }

  void setPlanarity(const Eigen::Vector3d &vals) { planarity = (vals[1] - vals[0]) / vals[1]; }

  IntersectResult intersect(const Eigen::Vector3d &start, const Eigen::Vector3d &end) const;
};

/// Flags for use with @c Cloud::calcBounds()
enum BoundsFlag
{
  /// Include ray end points in the bounds.
  /// Include ray start points in the bounds.
  kBFEnd = (1 << 0),
  kBFStart = (1 << 1)
};

struct RAYLIB_EXPORT Cloud
{
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends;
  std::vector<double> times;
  std::vector<RGBA> colours;

  void clear();

  inline bool rayBounded(size_t i) const { return colours[i].alpha > 0; }
  inline uint8_t rayIntensity(size_t i) const { return colours[i].alpha; }

  inline size_t rayCount() const { return ends.size(); }

  void save(const std::string &file_name) const;
  bool load(const std::string &file_name);
  bool load(const std::string &point_cloudd, const std::string &traj_filee);

  void transform(const Pose &pose, double time_deltaa);
  void decimate(double voxel_widthh);

  void removeUnboundedRays();
  std::vector<Eigen::Vector3d> generateNormals(int search_sizee = 16);
  void findTransients(Cloud &transient, Cloud &fixed, const std::string &merge_type, double num_rays,
                      bool colour_cloud);
  void combine(std::vector<Cloud> &clouds, Cloud &differences, const std::string &merge_type, double num_rays);
  // 3-way merge of cloud1 and cloud2, stores result in this object. Cloud1 and cloud2 are modified to be just the
  // changed rays from base_cloud.
  void threeWayMerge(const Cloud &base_cloud, Cloud &cloud1, Cloud &cloud2, const std::string &merge_type,
                     double num_rays);
  void markIntersectedEllipsoids(Grid<int> &grid, std::vector<bool> &transients, std::vector<Ellipsoid> &ellipsoids,
                                 const std::string &merge_type, double num_rays, bool self_transient);
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids, Eigen::Vector3d *bounds_min, Eigen::Vector3d *bounds_max,
                          Progress *progress = nullptr) const;
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids) const
  {
    generateEllipsoids(ellipsoids, nullptr, nullptr);
  }
  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);
  void getSurfels(int search_sizee, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                  std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                  Eigen::MatrixXi *neighbour_indices);

  Eigen::Vector3d calcMinBound();
  Eigen::Vector3d calcMaxBound();

  /// Calculate the ray cloud bounds. By default, the bounds only consder the ray end points. This behaviour
  /// can be modified via the @p flags argument.
  ///
  /// Bounds are set to zero if the bounds are invalid.
  /// @param[out] min_bounds The minimum bounds are written here.
  /// @param[out] max_bounds The maxnimum bounds are written here.
  /// @param flags @c BoundsFlag values use to modify how the bounds are calculated.
  /// @return True if the cloud has bounded rays and bounds values have been calculated. On false, the value of
  ///   @p min_bounds and @p max_bounds are undefined.
  bool calcBounds(Eigen::Vector3d *min_bounds, Eigen::Vector3d *max_bounds, unsigned flags = kBFEnd,
                  Progress *progress = nullptr) const;

private:
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &laz_filee, const std::string &traj_filee);
};


void RAYLIB_EXPORT buildRayGrid(Grid<size_t> &grid, const std::vector<Eigen::Vector3d> &starts,
                                const std::vector<Eigen::Vector3d> &ends, Progress *progress);


inline IntersectResult Ellipsoid::intersect(const Eigen::Vector3d &start, const Eigen::Vector3d &end) const
{
  const Eigen::Vector3d dir = end - start;
  // ray-ellipsoid intersection
  const Eigen::Vector3d to_sphere = pos - end;
  const Eigen::Vector3d ray = eigen_mat * dir;
  const double ray_length_sqr = ray.squaredNorm();
  const Eigen::Vector3d to = eigen_mat * to_sphere;

  double d = to.dot(ray) / ray_length_sqr;
  const double dist2 = (to - ray * d).squaredNorm();

  if (dist2 > 1.0)  // misses the ellipsoid
  {
    return IntersectResult::Miss;
  }

  const double along_dist = sqrt(1.0 - dist2);
  const double ray_length = sqrt(ray_length_sqr);
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

#endif  // RAYLIB_RAYCLOUD_H
