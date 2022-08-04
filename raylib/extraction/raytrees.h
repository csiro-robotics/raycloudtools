// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYEXTRACT_TREES_H
#define RAYLIB_RAYEXTRACT_TREES_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"
#include "raysegment.h"

namespace ray
{
/// Converts an index in to a unique colour
/// only for ints up to 255*255*255-1   (as it leaves black as a special colour)
inline void convertIntToColour(int x, RGBA &colour)
{
  colour.red = colour.green = colour.blue = 0;
  x++;
  for (int i = 0; i < 24; ++i)
  {
    if (x & (1 << i))
    {
      int channel = i % 3;
      int offset = i / 3;
      uint8_t &ch = channel == 0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
      ch |= (uint8_t)(1 << (7 - offset));
    }
  }
}

/// Converts a colour to a unique integer (index) value
/// returns -1 for the special case of the colour black, otherwise the integer index that the colour represents
inline int convertColourToInt(const RGBA &colour)
{
  int result = 0;
  for (int i = 0; i < 24; ++i)
  {
    int channel = i % 3;
    int offset = i / 3;
    const uint8_t &ch = channel == 0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
    if (ch & (1 << (7 - offset)))
      result |= 1 << i;
  }
  return result - 1;
}

/// structure containing the parameters used in tree reconstruction
struct TreesParams
{
  TreesParams()
    : max_diameter(0.9)
    , min_diameter(0.02)
    , distance_limit(1.0)
    , height_min(2.0)
    , length_to_radius(140.0)
    , cylinder_length_to_width(4.0)
    , gap_ratio(2.5)
    , span_ratio(4.5)
    , gravity_factor(0.3)
    , radius_exponent(0.67)
    , linear_range(3.0)
    , grid_width(0.0)
    , segment_branches(false)
  {}
  double max_diameter;     // maximum tree diameter. Trees wider than this may be segmented into multiple trees
  double min_diameter;     // minimum branch diameter. Branches thinner than this are not reconstructed 
  double distance_limit;   // maximum distance between points that can count as connected
  double height_min;       // minimum height for a tree. Lower values are considered undergrowth and excluded
  double length_to_radius; // the taper gradient of branches
  double cylinder_length_to_width; // the slenderness of the branch segment cylinders
  double gap_ratio;        // points with a wider gap determine that a branch has become two
  double span_ratio;       // points that span a larger width determine that a branch has become two
  double gravity_factor;   // preferences branches that are less lateral, so penalises implausable horizontal branches
  double radius_exponent;  // default 0.67 see "Allometric patterns in Acer platanoides (Aceraceae) branches"
                           // in "Wind loads and competition for light sculpt trees into self-similar structures" they
                           // suggest a range from 0.54 up to 0.89
  double linear_range;     // number of metres that branch radius is linear
  double grid_width;       // used on a grid cell with overlap, to remove trees with a base in the overlap zone
  bool segment_branches;   // flag to output the ray cloud coloured by branch segment index rather than by tree index
};

/// The structure for a single (cylindrical) branch section
struct BranchSection
{
  BranchSection()
    : tip(0, 0, 0)
    , radius(0)
    , parent(-1)
    , id(-1)
    , max_distance_to_end(0.0)
  {}
  Eigen::Vector3d tip;
  double radius;
  int parent;
  int id;  // 0 based per tree
  double max_distance_to_end;
  std::vector<int> roots;  // root points
  std::vector<int> ends;
  std::vector<int> children;
};

/// The structure for a set of trees, stored as a list of (connected) branch sections
/// together with the function for their extrsction from a ray cloud
struct Trees
{
  /// Constructs the piecewise cylindrical tree structures from the input ray cloud @c cloud
  /// The ground @c mesh defines the ground and @params are used to control the reconstruction 
  Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose);

  /// save the trees representation to a text file
  bool save(const std::string &filename);

  /// The piecewise cylindrical represenation of all of the trees
  std::vector<BranchSection> sections;
};

}  // namespace ray
#endif  // RAYLIB_RAYEXTRACT_TREES_H