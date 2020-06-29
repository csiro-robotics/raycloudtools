// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREEGEN_H
#define RAYLIB_RAYTREEGEN_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"

namespace ray
{
#define pitchAngle \
  (30.0 * kPi / 180.0)  // alternative to above, defined by relative angle of the head of the two branches
#define splitAngle (45.0 * kPi / 180.0)  // angle in the Y shape, usually around 45 degrees
#define branchGradient 20.0              // length per radius

double RAYLIB_EXPORT getMainBranchAngle(double covariance_angle);
void RAYLIB_EXPORT fillBranchAngleLookup();

/// Random generation of semi-realistic trees. These are based on a self-similar branching structure
/// for the @c pitchAngle @c splitAngle and @c branchGradient constants, with an additional @c random_factor
struct RAYLIB_EXPORT TreeGen
{
  /// create the tree structure, and list of leaf points
  void make(const Eigen::Vector3d &root_pos, double trunk_radius, double random_factor = 0.0);  // 0 to 1
  /// create a set of rays covering the tree at a roughly uniform distribution
  void generateRays(double ray_density);

  /// the ray cloud attributes
  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }
  
  struct Branch
  {
    Eigen::Vector3d tip;
    double radius;
    int parent_index;
  };
  /// access the geometry of the tree as a list of branches
  const std::vector<Branch> &branches() const { return branches_; }
  /// access the leaves of the tree
  const std::vector<Eigen::Vector3d> leaves() const { return leaves_; }
  /// the position of the base of the tree trunk
  const Eigen::Vector3d &root() const { return root_; }
  
private:
  std::vector<Eigen::Vector3d> leaves_;
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
  std::vector<Branch> branches_;
  Eigen::Vector3d root_;

  void addBranch(int parent_index, Pose pose, double radius, double random_factor);
};

}  // namespace ray

#endif  // RAYLIB_RAYTREEGEN_H
