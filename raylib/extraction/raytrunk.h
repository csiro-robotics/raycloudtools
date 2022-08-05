// Copyright (c) 2021
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
// Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
static const double minimum_score = 40.0;
static const double trunk_height_to_width = 4.0;  // height extent relative to real diameter of branch
static const double boundary_radius_scale = 3.0;  // how much farther out is the expected boundary compared to real
                                                  // branch radius? Larger requires more space to declare it a branch

/// Structure defining a single trunk, as used in trunk extraction
struct Trunk
{
  Trunk();
  Eigen::Vector3d centre;
  double radius;
  double score;
  double last_score;
  double length;
  double actual_length;
  Eigen::Vector3d dir;
  int parent;
  double tree_score;
  double distance_to_ground;
  double ground_height;
  bool active;
  bool visited;

  /// fill in the overlapping @c points to the trunk using the @c grid of points
  void getOverlap(const Grid<Eigen::Vector3d> &grid, std::vector<Eigen::Vector3d> &points, double spacing);

  /// estimate the centre and direction of the trunk from the shape of the points
  void estimatePose(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk direction vector estimation from the points
  void updateDirection(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk centre estimation from the neighbouring points
  void updateCentre(const std::vector<Eigen::Vector3d> &points);

  /// Update the trunk radius estimation from the points, and calculate its score (goodness of fit)
  void updateRadiusAndScore(const std::vector<Eigen::Vector3d> &points);
};

}  // namespace ray
#endif  // RAYLIB_RAYBRANCH_H