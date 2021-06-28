// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforestgen.h"
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
        for (auto &tree : trees())
        {
          double d = (root - tree.root()).norm();
          if (d < (radius + tree.branches()[0].radius) * 10.0)
          {
            found = false;
            break;
          }
        }
      }
      TreeGen tree;
      trees().push_back(tree);
      trees().back().make(root, radius, params.random_factor);
    }
    rad /= 2.0;
    num_trees *= pow(2.0, params.dimension);
  }
}

bool ForestGen::makeFromFile(const std::string &filename, const ForestParams &params)
{
  std::ifstream ifs(filename.c_str(), std::ios::out);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << std::endl;
    return false;
  }  
  std::vector<Eigen::Vector3d> bases;
  std::vector<double> radii;
  bool trunks_only = false;
  while (!ifs.eof())
  {
    Eigen::Vector3d base;
    double radius;
    std::string line;
    std::getline(ifs, line);
    if (line.length() == 0 || line[0] == '#')
      continue;
    int num_commas = (int)std::count(line.begin(), line.end(), ',');
    if (num_commas == 3) // just the base
    {
      trunks_only = true;
      std::istringstream ss(line);
      for (int i = 0; i<4; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i<3)
          base[i] = std::stod(token.c_str());
        else
          radius = std::stod(token.c_str());
      }
      bases.push_back(base);
      radii.push_back(radius);
    }
    else
    {
      if (trunks_only)
      {
        std::cerr << "bad input, some rows are just trunks and others aren't. Will not process correctly." << std::endl;
        return false;
      }
      TreeGen tree;
      trees().push_back(tree);
      trees().back().makeFromString(line);
    }
  }
  if (trunks_only)
    make(bases, radii, params);
  return true;
}

void ForestGen::make(const std::vector<Eigen::Vector3d> &roots, const std::vector<double> &radii, 
                     const ForestParams &params)
{
  ASSERT(roots.size() == radii.size());
  for (size_t i = 0; i<roots.size(); i++)
  {
    TreeGen tree;
    trees().push_back(tree);
    trees().back().make(roots[i], radii[i], params.random_factor);
  }
}

void ForestGen::generateRays(double ray_density)
{
  for (auto &tree : trees()) 
    tree.generateRays(ray_density);
}

std::vector<Eigen::Vector3d> ForestGen::getCanopy()
{
  std::vector<Eigen::Vector3d> canopy;
  for (auto &tree : trees()) 
    canopy.insert(canopy.end(), tree.leaves().begin(), tree.leaves().end());

  return canopy;
}

std::vector<Eigen::Vector3d> ForestGen::getPointCloud()
{
  std::vector<Eigen::Vector3d> cloud;
  for (auto &tree : trees()) 
    cloud.insert(cloud.end(), tree.rayEnds().begin(), tree.rayEnds().end());

  return cloud;
}
} // ray