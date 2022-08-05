// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREEGEN_H
#define RAYLIB_RAYTREEGEN_H

#include "raylib/raylibconfig.h"
#include "raytreestructure.h"
#include "raypose.h"
#include "rayutils.h"

namespace ray
{
#define pitchAngle \
  (30.0 * kPi / 180.0)  // alternative to above, defined by relative angle of the head of the two branches
#define splitAngle (45.0 * kPi / 180.0)  // angle in the Y shape, usually around 45 degrees
#define branchGradient 20.0              // length per radius

double RAYLIB_EXPORT getMainBranchAngle(double covariance_angle);
void RAYLIB_EXPORT fillBranchAngleLookup();

// store the basic parameters that describe the tree
struct TreeParams
{
  TreeParams() : min_branch_radius(0.001),
                 random_factor(0.0) {}
  double min_branch_radius;  
  double random_factor; // 0 to 1 value
};

/// Random generation of semi-realistic trees. These are based on a self-similar branching structure
/// for the @c pitchAngle @c splitAngle and @c branchGradient constants, with an additional @c random_factor
class RAYLIB_EXPORT TreeGen : public TreeStructure
{
public:
  TreeGen(){}
  TreeGen(const TreeStructure &base_tree)
  {
    segments_ = base_tree.segments();
    attribute_names_ = base_tree.attributes();
  }
  /// create the tree structure, and list of leaf points
  void make(const TreeParams &params);

  /// generate a set of rays as though the tree has been observed by a viewer circling it
  void generateRays(double ray_density);

  /// the rays generated from generateRays
  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }

  /// access the leaves of the tree
  const std::vector<Eigen::Vector3d> leaves() const { return leaves_; }  
private:
  std::vector<Eigen::Vector3d> leaves_;
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;

  void addBranch(int parent_index, Pose pose, double radius, const TreeParams &params);
};

}  // namespace ray

#endif  // RAYLIB_RAYTREEGEN_H
