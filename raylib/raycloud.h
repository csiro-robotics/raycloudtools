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

/// Flags for use with @c Cloud::calcBounds()
enum BoundsFlag
{
  /// Include ray end points in the bounds.
  /// Include ray start points in the bounds.
  kBFEnd = (1 << 0),
  kBFStart = (1 << 1)
};

/// This is the principle structure for representing a ray cloud. 
/// Rays are stored as line segments ( @c starts[i] to @c ends[i] ) together with a @c time and @c colour
/// The colour's alpha channel is used to store intensity, and so alpha=0 represents an unbounded ray
/// The end point of unbounded rays are not considered as part of the observed geometry of the scene
class RAYLIB_EXPORT Cloud
{
public:
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends;
  std::vector<double> times;
  std::vector<RGBA> colours;

  void clear();
  /// is the ray at index @c i bounded. Unbounded rays are non-returns, typically due to exceeding lidar range. 
  inline bool rayBounded(size_t i) const { return colours[i].alpha > 0; }
  /// this reflects the intensity of return recorded by the lidar. It is optional and does not affect the raycloudtools
  /// functions. However, the intensity should always be zero on any unbounded rays, as this information is used.
  inline uint8_t rayIntensity(size_t i) const { return colours[i].alpha; }

  /// the number of rays
  inline size_t rayCount() const { return ends.size(); }

  void save(const std::string &file_name) const;
  bool load(const std::string &file_name);
  bool load(const std::string &point_cloud, const std::string &traj_file);

  Eigen::Vector3d calcMinBound() const;
  Eigen::Vector3d calcMaxBound() const;

  /// apply a Euclidean transform and time shift to the ray cloud
  void transform(const Pose &pose, double time_delta);
  /// spatial decimation of the ray cloud, into one end point per voxel of width @c voxel_width
  void decimate(double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> *voxel_set = NULL);
  /// add a new ray to the ray cloud
  void addRay(const Eigen::Vector3d &start, const Eigen::Vector3d &end, double time, const RGBA &colour);
  /// add a new ray to the ray cloud, from another cloud
  void addRay(const Cloud &other_cloud, size_t index);

  void removeUnboundedRays();

  /// generates a covariance matrix of the nearest end points around each ray end in the cloud. The pointer arguments
  /// are optional attributes of this covariance matrix, which can be returned. Each covariance matrix represents a 
  /// SURFace ELement (surfel) with a centroid, normal, matrix and dimensions (of the ellipsoid that it represents)
  /// The list of neighbours can also be returned, to allow further analysis.
  void getSurfels(int search_size, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                  std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                  Eigen::MatrixXi *neighbour_indices);
  
  /// generates just the normal vectors of the ray end points based on each point's nearest neighbours.
  std::vector<Eigen::Vector3d> generateNormals(int search_size = 16);

  /// split a cloud based on the passed in function
  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);

  /// estimate the average spacing between end points of the ray cloud. This should be similar to the voxel
  /// width used on any spatially decimated ray clouds
  double estimatePointSpacing() const;

  /// Calculate the ray cloud bounds. By default, the bounds only consder the ray end points. This behaviour
  /// can be modified via the @p flags argument.
  ///
  /// Bounds are set to zero if the bounds are invalid.
  /// @param[out] min_bounds The minimum bounds are written here.
  /// @param[out] max_bounds The maximum bounds are written here.
  /// @param flags @c BoundsFlag values use to modify how the bounds are calculated.
  /// @return True if the cloud has bounded rays and bounds values have been calculated. On false, the value of
  ///   @p min_bounds and @p max_bounds are undefined.
  bool calcBounds(Eigen::Vector3d *min_bounds, Eigen::Vector3d *max_bounds, unsigned flags = kBFEnd,
                  Progress *progress = nullptr) const;

private:
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUD_H