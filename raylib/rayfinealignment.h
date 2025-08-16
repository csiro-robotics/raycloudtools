// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFINELIGNMENT_H
#define RAYLIB_RAYFINELIGNMENT_H

#include "raylib/raylibconfig.h"

#include "raycloud.h"
#include "rayutils.h"


namespace ray
{
/// Class for fine alignment of two ray clouds.
/// Being a gradient-descent based method, it requires the ray clouds to be nearly aligned at the start.
/// This means that the nearest surface on one ray cloud should be corresponding surface on the other cloud most of the
/// time.
class RAYLIB_EXPORT FineAlignment
{
public:
  /// Constructor takes two clouds as input @c clouds, also:
  /// @c non_rigid denotes whether the alignment transformation is quadratic or linear (Euclidean)
  /// @c verbose outputs debug text
  FineAlignment(Cloud *clouds, bool non_rigid, bool verbose)
    : clouds_(clouds)
    , non_rigid_(non_rigid)
    , verbose_(verbose)
  {}

  /// This function modifies clouds[0] (supplied in constructor) to match clouds[1]
  /// The alignment is either a rigid (Euclidean) transformation, or it contains some quadratic components to account
  /// for slight bend or warping within the cloud.
  void align();

private:
  /// Surfel object, suited to this alignment method
  struct RAYLIB_EXPORT Surfel
  {
    Surfel() {}
    Surfel(const Eigen::Vector3d &centroid, const Eigen::Matrix3d &matrix, const Eigen::Vector3d &width,
           const Eigen::Vector3d &normal, bool is_plane)
      : centroid(centroid)
      , matrix(matrix)
      , width(width)
      , normal(normal)
      , is_plane(is_plane)
    {}

    Eigen::Vector3d centroid;
    Eigen::Matrix3d matrix;
    Eigen::Vector3d width;
    Eigen::Vector3d normal;
    bool is_plane;
  };

  /// Identify matches between surfels by ID
  struct Match
  {
    int ids[2];
    Eigen::Vector3d normal;
  };

  /// A simple linear system structure. For solving Ax=b in least squares form (as AtA=Atb where t is transposition).
  struct LinearSystem
  {
    static const int state_size = 12;
    LinearSystem()
    {
      At_A.setZero();
      At_b.setZero();
    }
    Eigen::Matrix<double, state_size, 1> solve(bool verbose);
    Eigen::Matrix<double, state_size, state_size> At_A;
    Eigen::Matrix<double, state_size, 1> At_b;
  };

  /// Structure to store the nonlinear transformation
  struct QuadraticTransformation
  {
    QuadraticTransformation() {}
    QuadraticTransformation(const Eigen::Matrix<double, LinearSystem::state_size, 1> &x)
    {
      translation = Eigen::Vector3d(x[0], x[1], x[2]);
      rotation = Eigen::Vector3d(x[3], x[4], x[5]);
      a = Eigen::Vector3d(x[6], x[7], 0);
      b = Eigen::Vector3d(x[8], x[9], 0);
      c = Eigen::Vector3d(x[10], x[11], 0);
    }
    Eigen::Vector3d translation, rotation;  // the linear components
    Eigen::Vector3d a, b, c;                // the principle quadratic components
    /// This converts the linear translation and rotation vectors into a Pose (vector and quaternion) transformation.
    Pose getEuclideanPart() const
    {
      return Pose(translation, Eigen::Quaterniond(Eigen::AngleAxisd(rotation.norm(), rotation.normalized())));
    }
  };

  /// Create surfels per voxel of a vexelisation of the ray end points
  void generateSurfels();
  /// Find the list of correspondences between the two surfel sets surfels_[0] and surfels_[1]
  void generateSurfelMatches(std::vector<Match> &matches);
  /// Convert the matches into a linear system
  void buildLinearSystem(const std::vector<Match> &matches, double d, FineAlignment::LinearSystem &system);
  /// adjust the ray cloud 0 (and surfels_[0]) from the specified transformation @c trans
  void updateLinearSystem(std::vector<Match> &matches, const QuadraticTransformation &trans);

  /// Primary data:
  Cloud *clouds_;
  double non_rigid_;
  double verbose_;
  const double max_normal_difference_ = 0.5;
  bool no_normals_ {false};

  /// Derived data
  std::vector<Surfel> surfels_[2];
  double translation_weight_;
  Eigen::Vector3d centres_[2];
};
}  // namespace ray

#endif  // RAYLIB_RAYFINELIGNMENT_H
