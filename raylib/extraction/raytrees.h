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
  Trees(Cloud &cloud, const Mesh &mesh, bool verbose);
  bool save(const std::string &filename);
  std::vector<BranchSection> sections;
};

} // namespace ray
#endif // RAYLIB_RAYEXTRACT_TREES_H