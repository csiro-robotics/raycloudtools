// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#ifndef RAYINDEXINGOCTREENODE_H
#define RAYINDEXINGOCTREENODE_H

#include "raylib/raylibconfig.h"

#include "rayoctreenode.h"

#include <vector>

namespace ray
{
class RAYLIB_EXPORT IndexingOctreeNode : public OctreeNode
{
public:
  IndexingOctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max, bool is_root);
  IndexingOctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,  //
                     OctreeNode *child = nullptr, int child_index = 0);
  ~IndexingOctreeNode();

  void addDatum(size_t datum);
  bool removeDatum(size_t datum);

  const std::vector<size_t> &indices() const { return indices_; }

  void swapIndices(std::vector<size_t> &other) { std::swap(indices_, other); }

protected:
  OctreeNode *createNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max, OctreeNode *child,
                         int child_index) override;

private:
  std::vector<size_t> indices_;
};
}  // namespace ray

#endif  // RAYINDEXINGOCTREENODE_H
