// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFORESTSTRUCTURE_H
#define RAYLIB_RAYFORESTSTRUCTURE_H

#include "raylib/raylibconfig.h"
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raytreestructure.h"
#include "rayutils.h"

namespace ray
{
/// Underlying structure for sets of trees
struct RAYLIB_EXPORT ForestStructure
{
  std::vector<TreeStructure> trees;

  bool load(const std::string &filename);
  bool save(const std::string &filename);
  bool trunksOnly() { return trees.size() > 0 && trees[0].segments().size() == 1; }
  Eigen::Array<double, 9, 1> getMoments() const;
  void splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside);
  void generateSmoothMesh(Mesh &mesh, int red_id, double red_scale,
                          double green_scale, double blue_scale, bool add_uvs = false);
  /// reindex the segments to remove any disconnected segments, and order from root to tips
  void reindex();
  std::vector<std::string> comments; // just header comments
};
}  // namespace ray
#endif  // RAYLIB_RAYFORESTSTRUCTURE_H
