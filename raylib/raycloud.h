// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUD_H
#define RAYLIB_RAYCLOUD_H

#include "raylib/raycuboid.h"
#include "raylib/raylibconfig.h"

#include <set>
#include "raygrid.h"
#include "raypose.h"
#include "rayutils.h"

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
  /// reserve the cloud's vectors
  void reserve(size_t size);
  /// resize the cloud's vectors
  void resize(size_t size);

  /// is the ray at index @c i bounded. Unbounded rays are non-returns, typically due to exceeding lidar range.
  inline bool rayBounded(size_t i) const { return colours[i].alpha > 0; }
  /// this reflects the intensity of return recorded by the lidar. It is optional and does not affect the raycloudtools
  /// functions. However, the intensity should always be zero on any unbounded rays, as this information is used.
  inline uint8_t rayIntensity(size_t i) const { return colours[i].alpha; }

  /// the number of rays
  inline size_t rayCount() const { return ends.size(); }

  void save(const std::string &file_name) const;
  /// load a ray cloud file. @c check_extension checks the file extension before proceeding
  bool load(const std::string &file_name, bool check_extension = true, int min_num_rays = 4);

  /// minimum bounds of all bounded rays
  Eigen::Vector3d calcMinBound() const;
  /// maximum bounds of all bounded rays
  Eigen::Vector3d calcMaxBound() const;

  Eigen::Vector3d removeStartPos() // to aid in floating point accuracy 
  {
    Eigen::Vector3d offset(0,0,0);
    if (!ends.empty())
    {
      offset = ends[0];
      for (size_t i = 0; i<ends.size(); i++)
      {
        ends[i] -= offset;
        starts[i] -= offset;
      }
    }
    return offset;
  }
  void translate(const Eigen::Vector3d &offset)
  {
    if (offset.squaredNorm() != 0.0)
    {
      for (size_t i = 0; i<ends.size(); i++)
      {
        ends[i] += offset;
        starts[i] += offset;
      }
    }
  }

  /// apply a Euclidean transform and time shift to the ray cloud
  void transform(const Pose &pose, double time_delta);
  /// spatial decimation of the ray cloud, into one end point per voxel of width @c voxel_width
  void decimate(double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> &voxel_set);
  /// add a new ray to the ray cloud
  void addRay(const Eigen::Vector3d &start, const Eigen::Vector3d &end, double time, const RGBA &colour);
  /// add a new ray to the ray cloud, from another cloud
  void addRay(const Cloud &other_cloud, size_t index);

  void removeUnboundedRays();

  /// generates a covariance matrix of the nearest end points around each ray end in the cloud. The pointer arguments
  /// are optional attributes of this covariance matrix, which can be returned. Each covariance matrix represents a
  /// SURFace ELement (surfel) with a centroid, normal, matrix and dimensions (of the ellipsoid that it represents)
  /// The list of neighbours can also be returned, to allow further analysis.
  /// The last argument excludes back-facing rays from the surfel, this produces flatter surfels on thin double walls
  void getSurfels(int search_size, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                  std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                  Eigen::MatrixXi *neighbour_indices, double max_distance = 0.0,
                  bool reject_back_facing_rays = true) const;
  /// Get first and second order moments of cloud. This can be used as a simple way to compare clouds
  /// numerically. Note that different stats guarantee different clouds, but same stats do not guarantee same clouds
  /// These stats are arranged as: start mean, start sigma, end mean, end sigma, colour mean, time mean, time sigma,
  /// colour sigma. Times are doubles and colours vector4s, the others are vector3s, to give a total of 22 real values.
  Eigen::Array<double, 22, 1> getMoments() const;

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

  /// Static functions. These operate on the cloud file, and so do not require the full file to fit in memory

  /// Version for estimating the spacing between points for raycloud files.
  static double estimatePointSpacing(const std::string &file_name, const Cuboid &bounds, int num_points);

  /// Calculate the key information of a ray cloud, such as its bounds
  /// @c ends are only the bounded ones. @c starts are for all rays
  /// @c rays is all rays, so using the minimum known length for unbounded rays
  struct Info
  {
    // Axis-aligned bounding boxes
    Cuboid ends_bound;    // just the end points (not including for unbounded rays)
    Cuboid starts_bound;  // all start points
    Cuboid rays_bound;    // all ray extents

    int num_bounded;
    int num_rays;
    double min_time;
    double max_time;
    Eigen::Vector3d centroid;
    Eigen::Vector3d start_pos, end_pos;
  };
  static bool RAYLIB_EXPORT getInfo(const std::string &file_name, Info &info);

  /// Reads a ray cloud from file, and calls the function for each ray
  /// This forwards the call to a function appropriate to the ray cloud file format
  static bool read(const std::string &file_name,
                   std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                      std::vector<double> &times, std::vector<RGBA> &colours)>
                     apply);

private:
  bool loadPLY(const std::string &file, int min_num_rays);
  // Convert the set of neighbouring indices into a eigen solution, which is an ellipsoid of best fit.
  inline void eigenSolve(const std::vector<int> &ray_ids, const Eigen::MatrixXi &indices, int index, int num_neighbours,
                         Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> &solver, Eigen::Vector3d &centroid) const;
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUD_H
