// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBRANCH_H
#define RAYLIB_RAYBRANCH_H

#include "../raycloud.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"

namespace ray
{
static const double minimum_score = 40.0;         // unitless value representing trunk length per mean radius error of points
static const double trunk_height_to_width = 4.0;  // height to diameter ratio of the cylinder to fit to the trunk
static const double boundary_radius_scale = 3.0;  // how much farther out is the expected boundary compared to real
                                                  // branch radius? Larger requires more space to declare it a branch

/// Structure defining a single trunk, as used by raytrunks in trunk extraction
struct RAYLIB_EXPORT Trunk
{
  Trunk();
  Eigen::Vector3d centre;
  double radius;
  double score;
  double last_score;
  double length;
  double actual_length;
  double ground_height;
  Eigen::Vector3d dir;
  int parent;
  bool active;

  /// return the overlapping points to the trunk using the @c grid of points
  std::vector<Eigen::Vector3d> getOverlappingPoints(const Grid<Eigen::Vector3d> &grid, double spacing);

  /// estimate the centre and direction of the trunk from the shape of the points
  void estimatePose(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk direction vector (dir vector along cylinder axis) estimation from the points
  void updateDirection(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk centre estimation from the neighbouring points
  void updateCentre(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk radius estimation from the points
  void updateRadius(const std::vector<Eigen::Vector3d> &points);

  /// update the trunk candidate's score (goodness of fit)
  void updateScore(const std::vector<Eigen::Vector3d> &points);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace ray
#endif  // RAYLIB_RAYBRANCH_H