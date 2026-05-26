// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tom Lowe
#ifndef RAYDIFFER_H
#define RAYDIFFER_H

#include "raylib/raylibconfig.h"
#include "raylib/rayutils.h"

namespace ray
{

/// Get vector of closest distances between the point pairs
/// And point marked with float_max in their x component are excluded
void RAYLIB_EXPORT getDistancesBetweenPoints(const std::vector<Eigen::Vector3f> &points1, const std::vector<Eigen::Vector3f> &points2, std::vector<float> &dists_to_cloud1, std::vector<float> &dists_to_cloud2);

/// Print statistics on the distances, and adjust the distance_threshold if it is 0 to be at the statistical shoulder (bi-uniform fit)
/// returns the similarity (percentage under the distance threshold)
double RAYLIB_EXPORT printDistanceStatistics(const std::vector<float> &dists_to_cloud1, const std::vector<float> &dists_to_cloud2, double &distance_threshold);

}  // namespace ray

#endif  // RAYDIFFER_H
