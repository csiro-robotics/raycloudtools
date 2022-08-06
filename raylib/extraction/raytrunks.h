// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBRANCHES_H
#define RAYLIB_RAYBRANCHES_H

#include <iostream>
#include <map>
#include "../raycloud.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"
#include "raytrunk.h"

namespace ray
{
/// Structure storing a list of tree trunk descriptors. The constructor estimates
/// the trunks from the ray cloud
/// This is the class used by the tool rayextract trunks, which converts a ray cloud into
/// a text file describing its observed tree trunks
struct Trunks
{
  /// Reconstruct the set of trunks from the input ray cloud @c cloud, given a mean
  /// trunk radius @c midRadius.
  Trunks(const Cloud &cloud, double midRadius, bool verbose, bool exclude_passing_rays);

  /// Save the trunks to a text file
  bool save(const std::string &filename);

  /// Load the trunks from a text file
  static std::vector<std::pair<Eigen::Vector3d, double>> load(const std::string &filename);

  std::vector<Trunk> best_trunks;
};

}  // namespace ray
#endif  // RAYLIB_RAYBRANCHES_H