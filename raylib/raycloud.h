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

class RAYLIB_EXPORT Cloud
{
public:
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

  Eigen::Vector3d calcMinBound() const;
  Eigen::Vector3d calcMaxBound() const;

  void transform(const Pose &pose, double time_deltaa);
  void decimate(double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> *voxel_set = NULL);
  void addRay(const Eigen::Vector3d &start, const Eigen::Vector3d &end, double time, const RGBA &colour);
  void addRay(const Cloud &other_cloud, size_t index);

  void removeUnboundedRays();

  void getSurfels(int search_size, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                  std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                  Eigen::MatrixXi *neighbour_indices);
  std::vector<Eigen::Vector3d> generateNormals(int search_size = 16);

  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);
  double estimatePointSpacing() const;

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

}  // namespace ray

#endif  // RAYLIB_RAYCLOUD_H
