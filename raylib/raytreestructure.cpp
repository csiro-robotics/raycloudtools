#include "raytreestructure.h"

namespace ray
{
// calculate the tree's volume
double TreeStructure::volume()
{
  double volume = 0.0;
  for (size_t i = 1; i<segments_.size(); i++)
  {
    auto &branch = segments_[i];
    // fairly simple cylinder volume calculation...
    volume += (branch.tip - segments_[branch.parent_id].tip).norm() * branch.radius*branch.radius;
  }
  return volume * kPi; // .. but we multiply by pi at the end
}
} // namespace ray