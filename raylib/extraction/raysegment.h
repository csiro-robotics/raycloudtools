// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYSEGMENT_H
#define RAYLIB_RAYSEGMENT_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"


namespace ray
{

/// Vertex structure used in shortest path algorithm
struct RAYLIB_EXPORT Vertex
{
  Vertex(const Eigen::Vector3d &pos, const Eigen::Vector3d &start, uint8_t weight)
    : pos(pos), start(start)
    , parent(-1)
    , root(-1)
    , distance_to_ground(std::numeric_limits<double>::max())
    , distance_to_end(0.0)
    , score(std::numeric_limits<double>::max())
    , visited(false)
    , weight(weight)
  {}
  Eigen::Vector3d pos;
  Eigen::Vector3d start;
  int parent, root;
  double distance_to_ground;
  double distance_to_end;
  double score;
  bool visited;
  uint8_t weight;
};

/// Converts a ray cloud to a set of points @c points connected by the shortest path to the ground @c mesh
/// the returned vector of index sets provides the root points for each separated tree
/// @c height_min minimum height that counts as a tree
/// @c max_diameter maximum diameter of a tree trunk
/// @c distance_limit maximum distance between points that can be connected
/// @c gravity_factor controls how far laterally the shortest paths can travel
std::vector<std::vector<int>> RAYLIB_EXPORT getRootsAndSegment(std::vector<Vertex> &points, const Cloud &cloud, const Mesh &mesh,
                                                               double max_diameter, double distance_limit, double height_min,
                                                               double gravity_factor, bool alpha_weighting);

}  // namespace ray
#endif  // RAYLIB_RAYSEGMENT_H