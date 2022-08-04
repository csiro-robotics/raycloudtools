// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREEGEN_H
#define RAYLIB_RAYTREEGEN_H

#include "raylib/raylibconfig.h"

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

struct TreeParams
{
  TreeParams() : min_branch_radius(0.001),
                 random_factor(0.0) {}
  double min_branch_radius;  
  double random_factor; // 0 to 1 value
};

/// Random generation of semi-realistic trees. These are based on a self-similar branching structure
/// for the @c pitchAngle @c splitAngle and @c branchGradient constants, with an additional @c random_factor
class RAYLIB_EXPORT TreeStructure
{
public:
  /// create the tree structure, and list of leaf points
  void make(const TreeParams &params);
  void generateRays(double ray_density);

  /// the ray cloud attributes
  inline const std::vector<Eigen::Vector3d> rayStarts() const { return ray_starts_; }
  inline const std::vector<Eigen::Vector3d> rayEnds() const { return ray_ends_; }
  
  struct Segment
  {
    Segment() : tip(0,0,0), radius(0), parent_id(-1) {}
    Eigen::Vector3d tip;
    double radius;
    int parent_id;
    std::vector<double> attributes;
  };

  /// access the geometry of the tree as a list of branches
  const std::vector<Segment> &segments() const { return segments_; }
  std::vector<Segment> &segments() { return segments_; }
  /// access the leaves of the tree
  const std::vector<Eigen::Vector3d> leaves() const { return leaves_; }
  /// the position of the base of the tree trunk
  const Eigen::Vector3d &root() const { return segments_[0].tip; }
  std::vector<std::string> &attributes(){ return attribute_names_; }
  const std::vector<std::string> &attributes() const { return attribute_names_; }
  double volume();
  double &radius(){ return segments_[0].radius; }
  const double &radius() const { return segments_[0].radius; }
  
private:
  std::vector<Eigen::Vector3d> leaves_;
  std::vector<Eigen::Vector3d> ray_starts_, ray_ends_;
  std::vector<Segment> segments_;
  std::vector<std::string> attribute_names_;

  void addBranch(int parent_index, Pose pose, double radius, const TreeParams &params);
};

}  // namespace ray

#endif  // RAYLIB_RAYTREEGEN_H
