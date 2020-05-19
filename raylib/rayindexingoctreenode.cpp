// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "rayindexingoctreenode.h"

using namespace ray;

IndexingOctreeNode::IndexingOctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,
                                       bool is_root)
  : OctreeNode(bounds_min, bounds_max, is_root)
{}

IndexingOctreeNode::IndexingOctreeNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,  //
                                       OctreeNode *child, int child_index)
  : OctreeNode(bounds_min, bounds_max, child, child_index)
{}

IndexingOctreeNode::~IndexingOctreeNode() = default;

void IndexingOctreeNode::addDatum(size_t datum)
{
  indices_.emplace_back(datum);
}


bool IndexingOctreeNode::removeDatum(size_t datum)
{
  for (size_t i = 0; i < indices_.size(); ++i)
  {
    if (indices_[i] == datum)
    {
      std::swap(indices_[i], indices_[indices_.size() - 1]);
      indices_.pop_back();
      return true;
    }
  }
  return false;
}

OctreeNode *IndexingOctreeNode::createNode(const Eigen::Vector3d &bounds_min, const Eigen::Vector3d &bounds_max,
                                           OctreeNode *child, int child_index)
{
  return new IndexingOctreeNode(bounds_min, bounds_max, child, child_index);
}
