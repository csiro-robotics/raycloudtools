// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "rayoctreenode.h"

// #pragma GCC optimize("O0")

#include "rayintersect.h"

using namespace ray;

namespace
{
bool needsExpansion(const OctreeNode &tree, const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max,
                    Eigen::Vector3d *target)
{
  // We can't expand towards the AABB centre as that may be inside the tree. Instead we just expand towards the min
  // first, then to the max.
  if (!tree.contains(aabb_min))
  {
    *target = aabb_min;
    return true;
  }
  if (!tree.contains(aabb_max))
  {
    *target = aabb_max;
    return true;
  }
  return false;
}
}  // namespace

OctreeNode::OctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max, OctreeNode *child,
                       int child_index)
  : bounds_min_(bounds_min)
  , bounds_max_(bounds_max)
  , flags_(0)
{
  for (int i = 0; i < 8; ++i)
  {
    children_[i] = nullptr;
  }

  if (child)
  {
    children_[child_index] = child;
    // Child cannot be a root.
    child->flags_ &= ~F_Root;
    flags_ |= F_Branch;
  }
}

OctreeNode::OctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max, bool is_root)
  : OctreeNode(bounds_min, bounds_max, nullptr, 0)
{
  if (is_root)
  {
    flags_ |= F_Root;
  }
}

OctreeNode::~OctreeNode()
{
  for (int i = 0; i < 8; ++i)
  {
    delete children_[i];
  }
}


bool OctreeNode::contains(const Eigen::Vector3d &point) const
{
  return intersect::pointInAabb(point, bounds_min_, bounds_max_);
}


bool OctreeNode::overlaps(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max) const
{
  return intersect::aabbOverlap(bounds_min_, bounds_max_, aabb_min, aabb_max);
}


bool OctreeNode::rayIntersects(const Ray &ray) const
{
  double hit_times[2];
  if (!intersect::rayAabb(ray.origin, ray.direction, bounds_min_, bounds_max_, hit_times))
  {
    return false;
  }

  // Infinite ray intersects the box. Bounds check the ray.
  for (int i = 0; i < 2; ++i)
  {
    if (0 <= hit_times[i] && hit_times[i] < ray.length)
    {
      return true;
    }
  }

  return false;
}


OctreeNode *OctreeNode::addAabb(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max, double min_voxel_size)
{
  if (!contains(aabb_min) || !contains(aabb_max))
  {
    // This node is not big enough.
    return nullptr;
  }

  // Check to see if we are allowed to inspect children.
  // Allowed if this already is a branch or if child voxels will be at least min_voxel_size in size.
  if (isBranch() || (0.5 * (bounds_max_ - bounds_min_)).squaredNorm() >= min_voxel_size)
  {
    // First try find a child to contain the aabb.
    Eigen::Vector3d child_min, child_max;
    for (int i = 0; i < 8; ++i)
    {
      if (children_[i])
      {
        child_min = children_[i]->boundsMin();
        child_max = children_[i]->boundsMax();
      }
      else
      {
        boundsForChild(i, bounds_min_, bounds_max_, &child_min, &child_max);
      }

      if (intersect::pointInAabb(aabb_min, child_min, child_max) &&
          intersect::pointInAabb(aabb_max, child_min, child_max))
      {
        // Found a child which contains this box. Recurs on that child.
        if (!children_[i])
        {
          children_[i] = createNode(child_min, child_max, nullptr, 0);
          // Mark as branch.
          flags_ |= F_Branch;
        }
        return children_[i]->addAabb(aabb_min, aabb_max, min_voxel_size);
      }
    }
  }

  // Add to this node.
  return this;
}


void OctreeNode::rayTrace(const Ray &ray, const RayVisitFunction &visit, unsigned flags) const
{
  for (int i = 0; i < 8; ++i)
  {
    if (children_[i] && children_[i]->rayIntersects(ray))
    {
      if ((flags & children_[i]->flags()) != 0)
      {
        // Flags indicate visit should be called.
        visit(ray, children_[i], this);
      }
      // Even if visit wasn't called, we still recurse.
      children_[i]->rayTrace(ray, visit, flags);
    }
  }
}

OctreeNode *OctreeNode::expand(const Eigen::Vector3d &aabb_min, const Eigen::Vector3d &aabb_max)
{
  if (!isRoot())
  {
    return nullptr;
  }

  OctreeNode *tree = this;
  Eigen::Vector3d target_point;
  while (needsExpansion(*tree, aabb_min, aabb_max, &target_point))
  {
    Eigen::Vector3d bounds_min = tree->boundsMin();
    Eigen::Vector3d bounds_max = tree->boundsMax();
    int child_index = 0;

    // Determine the expansion direction. This determines the child index of tree in the new parent.
    child_index = expandBounds(target_point, &bounds_min, &bounds_max);

    if (child_index < 0)
    {
      // Failed to assign to child.
      return nullptr;
    }

    // Create a new octree node containing tree.
    // The direction of expanstion is important.
    OctreeNode *child = tree;
    tree = createNode(bounds_min, bounds_max, child, child_index);
    // Mark the new root.
    tree->flags_ |= F_Root;
  }

  return tree;
}


void OctreeNode::split(const OnSplitFuncion &on_split)
{
  if (!isLeaf())
  {
    return;
  }

  // Calculate bounds for children.
  Eigen::Vector3d child_min, child_max;
  for (int i = 0; i < 8; ++i)
  {
    boundsForChild(i, bounds_min_, bounds_max_, &child_min, &child_max);
    children_[i] = createNode(child_min, child_max, nullptr, 0);
  }

  // Mark as branch.
  flags_ |= F_Branch;

  if (on_split)
  {
    on_split(this, children_);
  }
}


bool OctreeNode::boundsForChild(int child_index, const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,
                                Eigen::Vector3d *child_min, Eigen::Vector3d *child_max)
{
  const Eigen::Vector3d half_extents = 0.5 * (bounds_max - bounds_min);

  *child_min = bounds_min;
  *child_max = bounds_max;

  switch (child_index)
  {
  case 0:
    // Left, back, bottom
    child_max->x() = child_min->x() + half_extents.x();
    child_max->y() = child_min->y() + half_extents.y();
    child_max->z() = child_min->z() + half_extents.z();
    break;

  case 1:
    // Right, back, bottom
    child_min->x() = bounds_min.x() + half_extents.x();
    child_max->x() = bounds_max.x();
    child_max->y() = child_min->y() + half_extents.y();
    child_max->z() = child_min->z() + half_extents.z();
    break;

  case 2:
    // Left, front, bottom.
    child_max->x() = child_min->x() + half_extents.x();
    child_min->y() = bounds_min.y() + half_extents.y();
    child_max->y() = bounds_max.y();
    child_max->z() = child_min->z() + half_extents.z();
    break;

  case 3:
    // Right, front, bottom
    child_min->x() = bounds_min.x() + half_extents.x();
    child_max->x() = bounds_max.x();
    child_min->y() = bounds_min.y() + half_extents.y();
    child_max->y() = bounds_max.y();
    child_max->z() = child_min->z() + half_extents.z();
    break;

  case 4:
    // ---
    // Left, back, top
    child_max->x() = child_min->x() + half_extents.x();
    child_max->y() = child_min->y() + half_extents.y();
    child_min->z() = bounds_min.z() + half_extents.z();
    child_max->z() = bounds_max.z();
    break;

  case 5:
    // Right, back, top
    child_min->x() = bounds_min.x() + half_extents.x();
    child_max->x() = bounds_max.x();
    child_max->y() = child_min->y() + half_extents.y();
    child_min->z() = bounds_min.z() + half_extents.z();
    child_max->z() = bounds_max.z();
    break;

  case 6:
    // Left, front, top.
    child_max->x() = child_min->x() + half_extents.x();
    child_min->y() = bounds_min.y() + half_extents.y();
    child_max->y() = bounds_max.y();
    child_min->z() = bounds_min.z() + half_extents.z();
    child_max->z() = bounds_max.z();
    break;

  case 7:
    // Right, front, top
    child_min->x() = bounds_min.x() + half_extents.x();
    child_max->x() = bounds_max.x();
    child_min->y() = bounds_min.y() + half_extents.y();
    child_max->y() = bounds_max.y();
    child_min->z() = bounds_min.z() + half_extents.z();
    child_max->z() = bounds_max.z();
    break;

  default:
    return false;
  }
  return true;
}


int OctreeNode::expandBounds(const Eigen::Vector3d &targetPoint, Eigen::Vector3d *bounds_min,
                             Eigen::Vector3d *bounds_max)
{
  const bool below_x = (targetPoint.x() < bounds_min->x());
  const bool below_y = (targetPoint.y() < bounds_min->y());
  const bool below_z = (targetPoint.z() < bounds_min->z());

  // Child indexing is:
  // 0: -x, -y, -z
  // 1:  x, -y, -z
  // 2: -x,  y, -z
  // 3:  x,  y, -z
  // 4: -x, -y,  z
  // 5:  x, -y,  z
  // 6: -x,  y,  z
  // 7:  x,  y,  z

  const Eigen::Vector3d bounds_extents = *bounds_max - *bounds_min;

  if (below_x && below_y && below_z)
  {
    // Expand such that we place the current node into the right, front, top.
    bounds_min->x() -= bounds_extents.x();
    bounds_min->y() -= bounds_extents.y();
    bounds_min->z() -= bounds_extents.z();
    return 7;
  }
  else if (!below_x && below_y && below_z)
  {
    // Expand from left front top.
    bounds_max->x() += bounds_extents.x();
    bounds_min->y() -= bounds_extents.y();
    bounds_min->z() -= bounds_extents.z();
    return 6;
  }
  else if (below_x && !below_y && below_z)
  {
    // Expand from right back top.
    bounds_min->x() -= bounds_extents.x();
    bounds_max->y() += bounds_extents.y();
    bounds_min->z() -= bounds_extents.z();
    return 5;
  }
  else if (!below_x && !below_y && below_z)
  {
    // Expand from left back top.
    bounds_max->x() += bounds_extents.x();
    bounds_max->y() += bounds_extents.y();
    bounds_min->z() -= bounds_extents.z();
    return 4;
  }
  if (below_x && below_y && !below_z)
  {
    // Expand from right, front, top.
    bounds_min->x() -= bounds_extents.x();
    bounds_min->y() -= bounds_extents.y();
    bounds_max->z() += bounds_extents.z();
    return 3;
  }
  else if (!below_x && below_y && !below_z)
  {
    // Expand from left front top.
    bounds_max->x() += bounds_extents.x();
    bounds_min->y() -= bounds_extents.y();
    bounds_max->z() += bounds_extents.z();
    return 2;
  }
  else if (below_x && !below_y && !below_z)
  {
    // Expand from right back top.
    bounds_min->x() -= bounds_extents.x();
    bounds_max->y() += bounds_extents.y();
    bounds_max->z() += bounds_extents.z();
    return 1;
  }
  else if (!below_x && !below_y && !below_z)
  {
    // Expand from left back top.
    bounds_max->x() += bounds_extents.x();
    bounds_max->y() += bounds_extents.y();
    bounds_max->z() += bounds_extents.z();
    return 0;
  }

  return -1;
}
