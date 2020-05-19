// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYOCTREENODE_H
#define RAYOCTREENODE_H

#include "raylib/raylibconfig.h"

#include <Eigen/Dense>

#include <cinttypes>
#include <functional>

namespace ray
{
class RAYLIB_EXPORT OctreeNode
{
public:
  struct RAYLIB_EXPORT Ray
  {
    Eigen::Vector3d origin;
    Eigen::Vector3d direction;
    Eigen::Vector3d end;
    double length;
    double length_inverse;

    inline Ray() {}
    Ray(const Eigen::Vector3d &origin, const Eigen::Vector3d &end);
    Ray(const Eigen::Vector3d &origin, const Eigen::Vector3d &direction, double length);
  };

  enum Flag : unsigned
  {
    F_Zero = 0,
    /// Flag marking the root node.
    F_Root = (1 << 0),
    /// Flag marking a ranch node. May also be the root.
    F_Branch = (1 << 1),
    /// Set if the node has data attached to it.
    F_HasData = (1 << 2),
    F_All = 0xffffu
  };

  /// Function prototype called when visiting a node in a ray trace.
  /// @param node The parent node being visited.
  /// @param parent Parent of @p node.
  using RayVisitFunction = std::function<void(const Ray &, const OctreeNode *, const OctreeNode *)>;
  using OnSplitFuncion = std::function<void(OctreeNode *, OctreeNode *const *)>;

  OctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max, bool is_root);
  OctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,  //
             OctreeNode *child = nullptr, int child_index = 0);
  virtual ~OctreeNode();

  inline bool isBranch() const { return flags_ & F_Branch; }
  inline bool isLeaf() const { return !isBranch(); }
  inline bool isRoot() const { return flags_ & F_Root; }
  inline bool hasData() const { return flags_ & F_HasData; }
  inline unsigned flags() const { return flags_; }

  inline Eigen::Vector3d boundsMin() const { return bounds_min_; }
  inline Eigen::Vector3d boundsMax() const { return bounds_max_; }

  inline OctreeNode *child(int index) { return children_[index]; }
  inline const OctreeNode *child(int index) const { return children_[index]; }

  bool contains(const Eigen::Vector3d &point) const;
  bool overlaps(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max) const;
  bool rayIntersects(const Ray &ray) const;

  OctreeNode *addAabb(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max, double min_voxel_size);

  /// Trace a ray recursively through this node calling @p visit for each intersected decendent node filtered by
  /// @p flags.
  ///
  /// The visit function is called only for descendents which have at least one flags matching @p flags.
  ///
  /// The Ray is assumed to intersect this node.
  void rayTrace(const Ray &ray, const RayVisitFunction &visit, unsigned flags = F_All) const;

  /// Creates a new OctreeNode which encompass this node and the given axis aligned bounding box (AABB).
  ///
  /// This method fails if @c isRoot() is false for this node.
  OctreeNode *expand(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max);

  /// Creates a new OctreeNode which encompass this node and the given @p target point.
  ///
  /// This method fails if @c isRoot() is false for this node.
  inline OctreeNode *expand(const Eigen::Vector3d &target) { return expand(target, target); }

  /// Split this (leaf) node into it's 8 children.
  ///
  /// After splitting, calls @p on_split to allow derivations to migrade data into the new leaf nodes.
  void split(const OnSplitFuncion &on_split = OnSplitFuncion());

  static bool boundsForChild(int child_index, const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,
                             Eigen::Vector3d *child_min, Eigen::Vector3d *child_max);

  static int expandBounds(const Eigen::Vector3d &targetPoint, Eigen::Vector3d *bounds_min, Eigen::Vector3d *bounds_max);

protected:
  /// Implemented by derivations to create new octree nodes.
  virtual OctreeNode *createNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,
                                 OctreeNode *child, int child_index) = 0;

private:
  Eigen::Vector3d bounds_min_;
  Eigen::Vector3d bounds_max_;
  OctreeNode *children_[8];
  unsigned flags_ = 0;
};


inline OctreeNode::Ray::Ray(const Eigen::Vector3d &origin, const Eigen::Vector3d &end)
  : origin(origin)
  , end(end)
{
  direction = end - origin;
  length = direction.norm();
  if (length > 1e-9)
  {
    length_inverse = 1.0 / length;
  }
  else
  {
    length_inverse = 0;
  }

  direction *= length_inverse;
}


inline OctreeNode::Ray::Ray(const Eigen::Vector3d &origin, const Eigen::Vector3d &direction, double length)
  : origin(origin)
  , direction(direction)
  , length(length)
{
  if (length > 1e-9)
  {
    length_inverse = 1.0 / length;
  }
  else
  {
    length_inverse = 0;
  }

  end = origin + direction * length;
}

}  // namespace ray

#endif  // RAYOCTREENODE_H
