// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREESUMMARY_H
#define RAYLIB_RAYTREESUMMARY_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"

namespace ray
{

struct TreeSummary
{
  Eigen::Vector3d base; 
  double height;
  double radius;
  bool trunk_identified;

  static std::vector<TreeSummary> load(const std::string &filename)
  {
    std::vector<TreeSummary> summaries;
    std::ifstream ifs(filename.c_str(), std::ios::out);
    if (!ifs.is_open())
    {
      std::cerr << "Error: cannot open " << filename << std::endl;
      return summaries;
    }  
    while (!ifs.eof())
    {
      std::string line;
      std::getline(ifs, line);
      if (line.length() == 0 || line[0] == '#')
        continue;
      int num_commas = (int)std::count(line.begin(), line.end(), ',');
      if (num_commas >= 3 && num_commas <= 5) 
      {
        std::istringstream ss(line);
        TreeSummary sum;
        sum.trunk_identified = true;
        for (int i = 0; i<num_commas+1; i++)
        {
          std::string token;
          std::getline(ss, token, ',');
          if (i<3)
            sum.base[i] = std::stod(token.c_str());
          else if (i==3)
            sum.radius = std::stod(token.c_str());
          else if (i==4)
            sum.height = std::stod(token.c_str());
          else if (i==5)
            sum.trunk_identified = (bool)std::stoi(token.c_str());
        }
        const double tree_height_per_radius = 20.0;
        if (num_commas == 3) // no height specified to estimate it
          sum.height = sum.radius * tree_height_per_radius;
        summaries.push_back(sum);
      }   
      else
      {
        std::cerr << "Unexpected file format, there should be 4 to 6 comma=separated variables" << std::endl;
        return summaries;
      } 
    }
    return summaries;
  }
  static bool save(const std::string &filename, const std::vector<TreeSummary> &summaries)
  {
    if (summaries.empty())
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
    ofs << "# Forest extraction, tree base location list: x, y, z, trunk radius, height, trunk_identified (radius is estimated from height when trunk not identified)" << std::endl;
    for (auto &result: summaries)
    {
      ofs << result.base[0] << ", " << result.base[1] << ", " << result.base[2] << ", " << result.radius << ", " << result.height << ", " << result.trunk_identified << std::endl;
    }
    return true;    
  }
};

}  // namespace ray

#endif  // RAYLIB_RAYTREESUMMARY_H
