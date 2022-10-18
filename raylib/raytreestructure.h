// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREESTRUCTURE_H
#define RAYLIB_RAYTREESTRUCTURE_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

namespace ray
{
/// Random generation of semi-realistic trees. These are based on a self-similar branching structure
/// for the @c pitchAngle @c splitAngle and @c branchGradient constants, with an additional @c random_factor
class RAYLIB_EXPORT TreeStructure
{
public:
  /// The tree is a list of segments, which are connected through the parent_id
  struct Segment
  {
    Segment()
      : tip(0, 0, 0)
      , radius(0)
      , parent_id(-1)
    {}
    Eigen::Vector3d tip;
    double radius;
    int parent_id;
    std::vector<double> attributes;
  };

  /// access the geometry of the tree as a list of branches
  const std::vector<Segment> &segments() const { return segments_; }
  std::vector<Segment> &segments() { return segments_; }
  /// the position of the base of the tree trunk
  const Eigen::Vector3d &root() const { return segments_[0].tip; }

  /// access the tree's attributes
  std::vector<std::string> &treeAttributes() { return tree_attribute_names_; }
  const std::vector<std::string> &treeAttributes() const { return tree_attribute_names_; }
  std::vector<std::string> &branchAttributes() { return branch_attribute_names_; }
  const std::vector<std::string> &branchAttributes() const { return branch_attribute_names_; }

  /// calculate the volume of the tree
  double volume();

  /// return the root radius of the tree
  double &radius() { return segments_[0].radius; }
  const double &radius() const { return segments_[0].radius; }

protected:
  std::vector<Segment> segments_;
  std::vector<std::string> tree_attribute_names_;
  std::vector<std::string> branch_attribute_names_;
};
}  // namespace ray
#endif  // RAYLIB_RAYTREESTRUCTURE_H
