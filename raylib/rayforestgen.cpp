// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforestgen.h"
#include "rayforeststructure.h"
#include "rayutils.h"

namespace ray
{
void ForestGen::make(const ForestParams &params)
{
  double rad = params.max_tree_radius;
  double num_trees = sqr(params.field_width) * params.adult_tree_density;
  for (int level = 0; level < 2; level++)
  {
    for (int i = 0; i < (int)num_trees; i++)
    {
      double radius = rad * (1.0 + random(-0.25, 0.5) * params.random_factor);
      Eigen::Vector3d root;
      bool found = false;
      while (!found)
      {
        root = params.field_width * 0.5 * Eigen::Vector3d(random(-1.0, 1.0), random(-1.0, 1.0), 0.0);
        found = true;
        for (auto &tree : trees_)
        {
          double d = (root - tree.root()).norm();
          if (d < (radius + tree.segments()[0].radius) * 10.0)
          {
            found = false;
            break;
          }
        }
      }
      TreeGen tree;
      TreeGen::Segment segment;
      segment.tip = root;
      segment.radius = radius;
      tree.segments().push_back(segment);
      trees_.push_back(tree);
      trees_.back().make(params);
    }
    rad /= 2.0;
    num_trees *= pow(2.0, params.dimension);
  }
}

bool ForestGen::makeFromFile(const std::string &filename, const TreeParams &params)
{
  ForestStructure forest;
  if (!forest.load(filename))
  {
    return false;
  }
  for (auto &tree : forest.trees)
  {
    trees_.push_back(TreeGen(tree));
  }
  if (trees_[0].segments().size() == 1)  // must have loaded trunks only
  {
    for (auto &tree : trees_)
    {
      tree.make(params);
    }
  }
  return true;
}

void ForestGen::generateRays(double ray_density)
{
  for (auto &tree : trees_)
  {
    tree.generateRays(ray_density);
  }
}

std::vector<Eigen::Vector3d> ForestGen::getCanopy()
{
  std::vector<Eigen::Vector3d> canopy;
  for (auto &tree : trees_)
  {
    canopy.insert(canopy.end(), tree.leaves().begin(), tree.leaves().end());
  }

  return canopy;
}

std::vector<Eigen::Vector3d> ForestGen::getPointCloud()
{
  std::vector<Eigen::Vector3d> cloud;
  for (auto &tree : trees_)
  {
    cloud.insert(cloud.end(), tree.rayEnds().begin(), tree.rayEnds().end());
  }

  return cloud;
}
}  // namespace ray