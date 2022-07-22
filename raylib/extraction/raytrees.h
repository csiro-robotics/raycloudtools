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
// only ints up to 255*255*255-1   (as it leaves black as a special colour)
inline void convertIntToColour(int x, RGBA &colour)
{
  colour.red = colour.green = colour.blue = 0;
  x++;
  for (int i = 0; i<24; ++i)
  {
    if (x & (1<<i))
    {
      int channel = i%3;
      int offset = i/3;
      uint8_t &ch = channel == 0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
      ch |= (uint8_t)(1<<(7-offset));
    }
  }
}

// returns -1 for the special case of the colour black, otherwise the integer index that the colour represents
inline int convertColourToInt(const RGBA &colour)
{
  int result = 0;
  for (int i = 0; i<24; ++i)
  {
    int channel = i%3;
    int offset = i/3;
    const uint8_t &ch = channel==0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
    if ( ch & (1<<(7-offset)))
      result |= 1<<i;
  }
  return result-1;
}


struct TreesParams
{
  TreesParams() : max_diameter(0.9), min_diameter(0.02), distance_limit(1.0), height_min(2.0), length_to_radius(140.0), 
                  cylinder_length_to_width(4.0), gap_ratio(2.5), span_ratio(4.5), gravity_factor(0.3), radius_exponent(0.67),
                  linear_range(3.0), grid_width(0.0), segment_branches(false) {}
  double max_diameter;
  double min_diameter;
  double distance_limit;
  double height_min;
  double length_to_radius;
  double cylinder_length_to_width;
  double gap_ratio;
  double span_ratio;
  double gravity_factor;
  double radius_exponent; // default 0.67 see "Allometric patterns in Acer platanoides (Aceraceae) branches"
                          // in "Wind loads and competition for light sculpt trees into self-similar structures" they
                          // suggest a range from 0.54 up to 0.89
  double linear_range; // number of metres that branch radius is linear
  double grid_width;
  bool segment_branches;
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