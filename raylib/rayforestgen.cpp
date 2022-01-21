// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforestgen.h"
#include "rayutils.h"

namespace ray
{
// Parse the tree file into a vector of tree structures
bool ForestStructure::load(const std::string &filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::out);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << std::endl;
    return false;
  }  
  std::string line;
  std::getline(ifs, line);
  std::string mandatory_text = "# tree file: x,y,z,radius";
  if (line.substr(0, mandatory_text.length()) != mandatory_text)
  {
    std::cerr << "Error: tree files must start with the mandatory files: " << mandatory_text << std::endl;
    return false;
  }
  int commas_per_segment = 1 + (int)std::count(line.begin(), line.end(), ',');
  std::istringstream ss(line);
  std::vector<std::string> attributes;
  bool has_parent_id = false;
  for (int i = 4; i<commas_per_segment; i++)
  {
    std::string token;
    std::getline(ss, token, ',');
    if (token == "parent_id")
      has_parent_id = true;
    else
      attributes.push_back(token);
  }

  int line_number = 0;
  while (!ifs.eof())
  {
    std::string line;
    std::getline(ifs, line);
    if (line.length() == 0 || line[0] == '#')
      continue;
    int num_commas = 1 + (int)std::count(line.begin(), line.end(), ',');
    TreeStructure tree;
    tree.attributes() = attributes;
    if (num_commas > commas_per_segment && !has_parent_id)
    {
      std::cerr << "Error: trees with multiple segments need parent_id field to connect them" << std::endl;
      return false;
    }
    if (num_commas%commas_per_segment != 0)
    {
      std::cerr << "Error: line " << line_number << " contains " << num_commas << " commas, which is not divisible by " << commas_per_segment << std::endl;
      return false;
    }
    line_number++;
    std::istringstream ss(line);
    for (int c = 0; c<num_commas; c += commas_per_segment)
    {
      TreeStructure::Segment segment; 
      for (int i = 0; i<commas_per_segment; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i<3)
          segment.tip[i] = std::stod(token.c_str());
        else if (i==3)
          segment.radius = std::stod(token.c_str());
        else if (i==4 && has_parent_id)
          segment.parent_id = std::stoi(token.c_str());
        else
          segment.attributes.push_back(std::stod(token.c_str()));
      }
      tree.segments().push_back(segment);
    }
    trees.push_back(tree);
  }
  return true;
}

bool ForestStructure::save(const std::string &filename)
{
  if (trees.empty())
  {
    std::cerr << "No data to save to " << filename << std::endl;
    return false;
  }
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree file: x,y,z,radius";
  if (trees[0].segments().size() > 1)
    ofs << ",parent_id";
  for (auto &att: trees[0].attributes())
    ofs << "," << att;
  ofs << std::endl;
  for (auto &tree: trees)
  {
    for (auto &segment: tree.segments())
    {
      ofs << segment.tip[0] << "," << segment.tip[1] << "," << segment.tip[2] << "," << segment.radius;
      if (tree.segments().size() > 1)
        ofs << "," << segment.parent_id;
      for (auto &att: segment.attributes)
        ofs << "," << att;
      ofs << " ";
    }
    ofs << std::endl;
  }
  return true;      
}


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
        for (auto &tree : trees)
        {
          double d = (root - tree.root()).norm();
          if (d < (radius + tree.segments()[0].radius) * 10.0)
          {
            found = false;
            break;
          }
        }
      }
      TreeStructure tree;
      TreeStructure::Segment segment;
      segment.tip = root;
      segment.radius = radius;
      tree.segments().push_back(segment);
      trees.push_back(tree);
      trees.back().make(params.random_factor);
    }
    rad /= 2.0;
    num_trees *= pow(2.0, params.dimension);
  }
}


bool ForestGen::makeFromFile(const std::string &filename, double random_factor)
{
  if (!load(filename))
    return false;
  for (auto &tree: trees)
    tree.make(random_factor);
  return true;
}

void ForestGen::generateRays(double ray_density)
{
  for (auto &tree : trees) 
    tree.generateRays(ray_density);
}

std::vector<Eigen::Vector3d> ForestGen::getCanopy()
{
  std::vector<Eigen::Vector3d> canopy;
  for (auto &tree : trees) 
    canopy.insert(canopy.end(), tree.leaves().begin(), tree.leaves().end());

  return canopy;
}

std::vector<Eigen::Vector3d> ForestGen::getPointCloud()
{
  std::vector<Eigen::Vector3d> cloud;
  for (auto &tree : trees) 
    cloud.insert(cloud.end(), tree.rayEnds().begin(), tree.rayEnds().end());

  return cloud;
}
} // ray