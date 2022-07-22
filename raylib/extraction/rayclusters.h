// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_CLUSTERS_H
#define RAYLIB_CLUSTERS_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"

namespace ray
{
std::vector< std::vector<int> > generateClusters(const std::vector<Eigen::Vector3d> &points, double min_diameter, double max_diameter, bool agglomerate, bool verbose = false);
} // namespace ray
#endif // RAYLIB_RAYEXTRACT_WOODS_H