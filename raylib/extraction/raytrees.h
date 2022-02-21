// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYEXTRACT_TREES_H
#define RAYLIB_RAYEXTRACT_TREES_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include "../raymesh.h"
#include "raysegment.h"

namespace ray
{
struct TreesParams
{
  TreesParams() : max_diameter(0.9), distance_limit(1.0), height_min(2.0), minimum_radius(0.02), length_to_radius(80.0), 
                  cylinder_length_to_width(4.0), gap_ratio(2.5), span_ratio(4.5), gravity_factor(0.3) {}
  double max_diameter;
  double distance_limit;
  double height_min;
  double minimum_radius;
  double length_to_radius;
  double cylinder_length_to_width;
  double gap_ratio;
  double span_ratio;
  double gravity_factor;
};

struct BranchSection
{
  BranchSection() : tip(0,0,0), radius(0), parent(-1), id(-1), max_distance_to_end(0.0) {}
  Eigen::Vector3d tip;
  double radius;
  int parent;
  int id; // 0 based per tree
  double max_distance_to_end;
  std::vector<int> roots; // root points
  std::vector<int> ends;
  std::vector<int> children;
};    

struct Trees
{
  Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose);
  bool save(const std::string &filename);
  std::vector<BranchSection> sections;
};

} // namespace ray
#endif // RAYLIB_RAYEXTRACT_TREES_H