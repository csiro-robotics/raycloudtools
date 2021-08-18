// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBRANCH_H
#define RAYLIB_RAYBRANCH_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"

namespace ray
{
struct Branch
{
  Branch();
  Eigen::Vector3d centre; 
  double radius;
  double score;
  double length; 
  Eigen::Vector3d dir;
  int parent;
  double tree_score;
  double distance_to_ground;
  bool active;
  bool visited;

  void getOverlap(const Grid<Eigen::Vector3d> &grid, std::vector<Eigen::Vector3d> &points, double spacing);
  void estimatePose(const std::vector<Eigen::Vector3d> &points);
  void updateDirection(const std::vector<Eigen::Vector3d> &points, bool trunks_only);
  void updateCentre(const std::vector<Eigen::Vector3d> &points);
  void updateRadiusAndScore(const std::vector<Eigen::Vector3d> &points, double spacing, bool trunks_only);
};
} // namespace ray
#endif // RAYLIB_RAYBRANCH_H