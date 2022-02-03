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
#include "raysegment.h"

namespace ray
{
struct BranchSection
{
  BranchSection() : tip(0,0,0), radius(0), parent(-1), id(-1), rank(0) {}
  Eigen::Vector3d tip;
  double radius;
  int parent;
  int id; // 0 based per tree
  int rank; // number of segments up from the bottom segment
  std::vector<int> roots; // root points
  std::vector<int> ends;
  std::vector<int> children;
};    

struct Trees
{
  Trees(const Cloud &cloud, const std::vector<std::pair<Eigen::Vector3d, double> > &trunks, bool verbose);
  bool save(const std::string &filename);
  std::vector<BranchSection> sections;
};

} // namespace ray
#endif // RAYLIB_RAYEXTRACT_TREES_H