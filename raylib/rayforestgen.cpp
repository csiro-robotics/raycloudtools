// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforestgen.h"
#include "rayutils.h"
// #define OUTPUT_MOMENTS

namespace ray
{
// Parse the tree file into a vector of tree structures
bool ForestStructure::load(const std::string &filename)
{
  std::cout << "loading tree file: " << filename << std::endl;
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << std::endl;
    return false;
  }  
  std::string line;
  do
  {
    std::getline(ifs, line);
  } while (line[0] == '#');
  std::string mandatory_text = "x,y,z,radius";
  if (line.substr(0, mandatory_text.length()) != mandatory_text)
  {
    std::cerr << "Error: tree files must start with the mandatory format: " << mandatory_text << std::endl;
    return false;
  }
  int commas_per_segment = 1 + (int)std::count(line.begin(), line.end(), ',');
  line = line.substr(mandatory_text.length());
  if (line[0] == ',')
    line = line.substr(1);
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
  if (attributes.size() > 0)
  {
    std::cout << "reading extra tree attributes: ";
    for (auto &at: attributes)
      std::cout << at << "; ";
    std::cout << std::endl;
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
    if (has_parent_id && tree.segments().size() == 1)
    {
      std::cerr << "Error: format contains parent id, but trunk only (single segment) detected on line " << line_number << std::endl;
      return false;
    }
    trees.push_back(tree);
  }
  if (trees.empty())
    return false;
  if (trees[0].segments().empty())
    return false;

#if defined OUTPUT_MOMENTS
  Eigen::Array<double, 6, 1> mom = getMoments();
  std::cout << "stats: " << std::endl;
  for (int i = 0; i<mom.rows(); i++)
    std::cout << ", " << mom[i];
  std::cout << std::endl;
#endif // defined OUTPUT_MOMENTS    
  return true;
}

bool ForestStructure::save(const std::string &filename)
{
  if (trees.empty())
  {
    std::cerr << "No data to save to " << filename << std::endl;
    return false;
  }
  std::cout << "outputting tree file: " << filename << std::endl;
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree file:" << std::endl;
  ofs << "x,y,z,radius";
  if (trees[0].segments().size() > 1)
    ofs << ",parent_id";
  if (trees[0].attributes().size() > 0)
    std::cout << "Saving additional attributes: ";
  for (auto &att: trees[0].attributes())
  {
    ofs << "," << att;
    std::cout << att << ", ";
  }
  std::cout << std::endl;
  ofs << std::endl;
  for (auto &tree: trees)
  {
    for (size_t i = 0; i<tree.segments().size(); i++)
    {
      auto &segment = tree.segments()[i];
      ofs << segment.tip[0] << "," << segment.tip[1] << "," << segment.tip[2] << "," << segment.radius;
      if (tree.segments().size() > 1)
        ofs << "," << segment.parent_id;
      for (auto &att: segment.attributes)
        ofs << "," << att;
      if (i != tree.segments().size()-1)
        ofs << ", ";
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
      trees.back().make(params);
    }
    rad /= 2.0;
    num_trees *= pow(2.0, params.dimension);
  }
}

Eigen::Array<double, 6, 1> ForestStructure::getMoments() const
{
  Eigen::Array<double, 6, 1> moments;
  moments[0] = (double)trees.size();
  Eigen::Vector3d sum(0,0,0);
  Eigen::Vector3d sum_sqr(0,0,0);
  double rad = 0;
  double rad_sqr = 0.0;
  double volume = 0.0;
  for (auto &tree: trees)
  {
    Eigen::Vector3d p = tree.segments()[0].tip;
    sum += p;
    sum_sqr += Eigen::Vector3d(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
    double r = tree.segments()[0].radius;
    rad += r;
    rad_sqr += r*r;
    for (auto &segment: tree.segments())
    {
      if (segment.parent_id != -1)
      {
        double length = (segment.tip - tree.segments()[segment.parent_id].tip).norm();
        double r = segment.radius;
        volume += length * r * r;
      }
    }
  }
  moments[1] = sum.norm();
  moments[2] = sum_sqr.norm();
  moments[3] = rad;
  moments[4] = rad_sqr;
  moments[5] = volume * kPi;
  return moments;
}


bool ForestGen::makeFromFile(const std::string &filename, const TreeParams &params)
{
  if (!load(filename))
    return false;
  if (trees[0].segments().size() == 1) // must have loaded trunks only
  {
    for (auto &tree: trees)
      tree.make(params);
  }
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
}  // namespace ray