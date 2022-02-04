// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYSEGMENT_H
#define RAYLIB_RAYSEGMENT_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include <queue>


namespace ray
{
struct QueueNode
{
  QueueNode(){}
  QueueNode(double distance_to_ground, double score, double radius, int root, int index) : distance_to_ground(distance_to_ground), score(score), radius(radius), root(root), id(index) {}

  double distance_to_ground;
  double score;
  double radius;
  int root;
  int id;
};

#define MINIMISE_SQUARE_DISTANCE // bad: end points are so distant that it creates separate branches
#define MINIMISE_ANGLE // works quite well in flowing along branches, but sometimes causes multi-branch problem, where radius was too small. 
class QueueNodeComparator 
{ 
public: 
    bool operator() (const QueueNode &p1, const QueueNode &p2) 
    { 
#if defined MINIMISE_SQUARE_DISTANCE || defined MINIMISE_ANGLE
        return p1.score > p2.score; 
#else
        return p1.distance_to_ground > p2.distance_to_ground; 
#endif
    } 
}; 

static constexpr double inf = 1e10;

struct Vertex
{
  Vertex(){}
  Vertex(const Eigen::Vector3d &pos) : pos(pos), edge_pos(0,0,0), parent(-1), root(-1), distance_to_ground(inf), score(inf), visited(false) {}
  Eigen::Vector3d pos;
  Eigen::Vector3d edge_pos;
  int parent, root;
  double distance_to_ground;
  double score;
  bool visited;
};

void connectPointsShortestPath(std::vector<Vertex> &points, std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node, double distance_limit);
void segmentTrees(Cloud &cloud, double max_diameter, double gradient, double distance_limit, double height_min);

} // namespace ray
#endif // RAYLIB_RAYSEGMENT_H