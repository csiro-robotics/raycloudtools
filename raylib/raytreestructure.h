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
  Eigen::Vector3d closestPointOnSegment(int index, const Eigen::Vector3d &pos, Eigen::Vector3d &line_closest)
  {
    auto &segment = segments_[index];
    if (segment.parent_id == -1) // then treat as spherical
    {
      return segment.tip + (pos - segment.tip).normalized() * segment.radius;
    }
    Eigen::Vector3d dif = segment.tip - segments_[segment.parent_id].tip;
    double d = (pos - segments_[segment.parent_id].tip).dot(dif)/dif.squaredNorm();
    d = std::max(0.0, std::min(d, 1.0));
    line_closest = segments_[segment.parent_id].tip + dif * d;
    return line_closest + (pos - line_closest).normalized() * segment.radius;
  }

  /// access the geometry of the tree as a list of branches
  const std::vector<Segment> &segments() const { return segments_; }
  std::vector<Segment> &segments() { return segments_; }
  /// access the whole tree's attributes
  const std::vector<double> &treeAttributes() const { return tree_attributes_; }
  std::vector<double> &treeAttributes() { return tree_attributes_; }

  /// the position of the base of the tree trunk
  const Eigen::Vector3d &root() const { return segments_[0].tip; }

  /// access the tree's attributes
  std::vector<std::string> &treeAttributeNames() { return tree_attribute_names_; }
  const std::vector<std::string> &treeAttributeNames() const { return tree_attribute_names_; }
  std::vector<std::string> &attributeNames() { return branch_attribute_names_; }
  const std::vector<std::string> &attributeNames() const { return branch_attribute_names_; }

  /// calculate the volume of the tree
  double volume() const;

  /// return the root radius of the tree
  double &radius() { return segments_[0].radius; }
  const double &radius() const { return segments_[0].radius; }
  /// reindex the segments to remove any disconnected segments, and order from root to tips
  void reindex();

protected:
  std::vector<double> tree_attributes_;
  std::vector<Segment> segments_;
  std::vector<std::string> tree_attribute_names_;
  std::vector<std::string> branch_attribute_names_;
};
}  // namespace ray
#endif  // RAYLIB_RAYTREESTRUCTURE_H
