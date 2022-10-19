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

int TreeStructure::numSharedAttributes() const
{
  int n = 0;
  for (; n<(int)treeAttributes().size(); n++)
  {
    if (n == (int)branchAttributes().size())
    {
      break;
    }
    if (treeAttributes()[n] != branchAttributes()[n])
    {
      break;
    }
  }
  return n;
}
}  // namespace ray