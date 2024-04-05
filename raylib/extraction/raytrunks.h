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
class RAYLIB_EXPORT Trunks
{
public:
  /// Reconstruct the set of trunks from the input ray cloud @c cloud, given a mean
  /// trunk radius @c midRadius.
  Trunks(const Cloud &cloud, const Eigen::Vector3d &offset, double midRadius, bool verbose, bool remove_permeable_trunks);

  /// Save the trunks to a text file
  bool save(const std::string &filename, const Eigen::Vector3d &offset) const;

  /// Load the trunks from a text file
  static std::vector<std::pair<Eigen::Vector3d, double>> load(const std::string &filename);

  /// render trunk points to disk:
  void saveDebugTrunks(const std::string &filename, bool verbose, const std::vector<int> &lowest_trunk_ids, 
    const std::vector<Trunk> &trunks, const Eigen::Vector3d &offset) const;

  /// a forest nearest path search to find only the lowest trunks of any connected chain
  std::vector<int> findLowestTrunks(const std::vector<Trunk> &trunks) const;

  /// set the ground heights for each trunk
  void setTrunkGroundHeights(const Cloud &cloud, std::vector<Trunk> &trunks, 
    const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound);

  /// remove trunk candidates with rays that pass right through them
  void removePermeableTrunks(bool verbose, const Cloud &cloud, const Eigen::Vector3d &offset, std::vector<Trunk> &trunks, 
    const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound);

private:
  std::vector<Trunk> best_trunks_;
};

}  // namespace ray
#endif  // RAYLIB_RAYBRANCHES_H