// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_CLUSTERS_H
#define RAYLIB_CLUSTERS_H

#include "../rayutils.h"
#include "raylib/raylibconfig.h"

namespace ray
{
/// Generate clusters from a set of points based on min and max diameter criteria
/// This is an agglomerative clustering method, whereby clusters iteratively grow and fuse
void generateClusters(std::vector<std::vector<int>> &point_clusters, const std::vector<Eigen::Vector3d> &points, 
  double min_diameter, double max_diameter, bool verbose = false);
}  // namespace ray
#endif  // RAYLIB_RAYEXTRACT_WOODS_H