// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytreestructure.h"

namespace ray
{
// calculate the tree's volume
double TreeStructure::volume() const
{
  double volume = 0.0;
  for (size_t i = 1; i < segments_.size(); i++)
  {
    auto &branch = segments_[i];
    // fairly simple cylinder volume calculation...
    volume += (branch.tip - segments_[branch.parent_id].tip).norm() * branch.radius * branch.radius;
  }
  return volume * kPi;  // .. but we multiply by pi at the end
}

// order the indices from root to leaf contiguously and remove any disconnected segments
void TreeStructure::reindex()
{
  // first get the upwards info
  std::vector<std::vector<int> > children(segments_.size());
  for (size_t i = 0; i<segments_.size(); i++)
  {
    if (segments_[i].parent_id != -1)
    {
      children[segments_[i].parent_id].push_back(static_cast<int>(i));
    }
  }
  // now start at the root and work upwards:
  std::vector<Segment> new_segments;
  new_segments.push_back(segments_[0]);
  std::vector<int> child_segs = {0};
  for (size_t i = 0; i<child_segs.size(); i++)
  {
    int id = child_segs[i];
    for (auto &child_id : children[id])
    {
      Segment child = segments_[child_id];
      child.parent_id = static_cast<int>(i);
      new_segments.push_back(child);
      child_segs.push_back(child_id);
    }
  }
  segments_ = new_segments;
}

}  // namespace ray