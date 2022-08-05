// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYSEGMENT_H
#define RAYLIB_RAYSEGMENT_H

#include <queue>
#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"


namespace ray
{
static constexpr double inf = 1e10;

/// Vertex structure used in shortest path algorithm
struct Vertex
{
  Vertex() {}
  Vertex(const Eigen::Vector3d &pos)
    : pos(pos)
    , edge_pos(0, 0, 0)
    , parent(-1)
    , root(-1)
    , distance_to_ground(inf)
    , distance_to_end(0.0)
    , score(inf)
    , visited(false)
  {}
  Eigen::Vector3d pos;
  Eigen::Vector3d edge_pos;
  int parent, root;
  double distance_to_ground;
  double distance_to_end;
  double score;
  bool visited;
};

/// Converts a ray cloud to a set of points @c points connected by the shortest path to the ground @c mesh
/// the returned vector of index sets provides the root points for each separated tree
/// @c height_min minimum height that counts as a tree
/// @c max_diameter maximum diameter of a tree trunk
/// @c distance_limit maximum distance between points that can be connected
/// @c gravity_factor controls how far laterally the shortest paths can travel
std::vector<std::vector<int>> getRootsAndSegment(std::vector<Vertex> &points, Cloud &cloud, const Mesh &mesh,
                                                 double max_diameter, double distance_limit, double height_min,
                                                 double gravity_factor);

}  // namespace ray
#endif  // RAYLIB_RAYSEGMENT_H