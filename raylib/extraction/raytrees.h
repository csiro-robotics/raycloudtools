// Copyright (c) 2022
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
/// structure containing the parameters used in tree reconstruction
struct RAYLIB_EXPORT TreesParams
{
  TreesParams();
  double max_diameter;              // maximum tree diameter. Trees wider than this may be segmented into multiple trees
  double min_diameter;              // minimum branch diameter. Branches thinner than this are not reconstructed
  double distance_limit;            // maximum distance between points that can count as connected
  double height_min;                // minimum height for a tree. Lower values are considered undergrowth and excluded
  double length_to_radius;          // the taper gradient of branches
  double cylinder_length_to_width;  // the slenderness of the branch segment cylinders
  double gap_ratio;                 // points with a wider gap determine that a branch has become two
  double span_ratio;                // points that span a larger width determine that a branch has become two
  double gravity_factor;   // preferences branches that are less lateral, so penalises implausable horizontal branches
  double radius_exponent;  // default 0.67 see "Allometric patterns in Acer platanoides (Aceraceae) branches"
                           // in "Wind loads and competition for light sculpt trees into self-similar structures" they
                           // suggest a range from 0.54 up to 0.89
  double linear_range;     // number of metres that branch radius is linear
  double grid_width;       // used on a grid cell with overlap, to remove trees with a base in the overlap zone
  bool segment_branches;   // flag to output the ray cloud coloured by branch segment index rather than by tree index
};

struct BranchSection;  // forwards declaration

/// The class for a set of trees, stored as a list of (connected) branch sections
/// together with the function for their extrsction from a ray cloud
/// This is the class used in the tool rayextrsact trees, which reconstructs the branch structures
/// of the specified ray cloud, and saves them as a text file
class RAYLIB_EXPORT Trees
{
public:
  /// Constructs the piecewise cylindrical tree structures from the input ray cloud @c cloud
  /// The ground @c mesh defines the ground and @params are used to control the reconstruction
  Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose);

  /// save the trees representation to a text file
  bool save(const std::string &filename) const;

private:
  /// The piecewise cylindrical represenation of all of the trees
  std::vector<BranchSection> sections_;

  /// estimate branch radius from its length
  double radFromLength(double length) const;
  /// calculate the distance to farthest connected branch tip, for each point in the cloud
  void calculatePointDistancesToEnd();
  /// create the start branch segments at the root positions
  void generateRootSections(const std::vector<std::vector<int>> &roots_list);
  /// finalise the attributes of an end (tip) of a branch
  void setBranchTip();
  /// get the root position for the current section
  Eigen::Vector3d getRootPosition() const;
  /// find the points and end points within this branch section
  void extractNodesAndEndsFromRoots(std::vector<int> &nodes, const Eigen::Vector3d &base,
                                    const std::vector<std::vector<int>> &children);
  /// find separate clusters of points within the branch section
  std::vector<std::vector<int>> findPointClusters(const Eigen::Vector3d &base, bool &points_removed);
  /// split the branch section to one branch for each cluster
  void bifurcate(const std::vector<std::vector<int>> &clusters);
  /// find the points within the branch section from its end points
  void extractNodesFromEnds(std::vector<int> &nodes);
  /// set the branch section tip position from the supplied list of Vertex IDs
  Eigen::Vector3d calculateTipFromVertices(const std::vector<int> &nodes) const;
  /// estimate the vector to the cylinder centre from the set of nodes
  Eigen::Vector3d vectorToCylinderCentre(const std::vector<int> &nodes, const Eigen::Vector3d &dir) const;
  /// estimate the cylinder's radius from its centre, @c dir and set of nodes
  double estimateCylinderRadius(const std::vector<int> &nodes, const Eigen::Vector3d &dir) const;
  /// add a new section to continue reconstructing the branch
  void addChildSection();
  /// calculate the ownership, what branch section does each point belong to
  void calculateSectionIds(const std::vector<std::vector<int>> &roots_list, std::vector<int> &section_ids,
                           const std::vector<std::vector<int>> &children);
  /// debug draw
  void drawTrees(bool verbose);
  /// set ids that are locel (0-based) per tree
  void generateLocalSectionIds();
  /// if using an overlapping grid, then remove trees with base outside the non-overlapping cell bounds
  void removeOutOfBoundSections(const Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound);
  /// colour the cloud based on the section id for each point
  void segmentCloud(Cloud &cloud, std::vector<int> &root_segs, const std::vector<int> &section_ids);
  /// remove points from the ray cloud if outside of the non-overlapping grid cell bounds
  void removeOutOfBoundRays(Cloud &cloud, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound,
                            const std::vector<int> &root_segs);
  // cached data that is used throughout the processing method
  size_t sec_;
  const TreesParams *params_;
  std::vector<Vertex> points_;
  double max_radius_;
  double radius_length_scale_; /// Cached ratio from branch taper parameters
};

/// The structure for a single (cylindrical) branch section
struct RAYLIB_EXPORT BranchSection
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

/// Converts an index in to a unique colour
/// only for ints up to 255*255*255-1   (as it leaves black as a special colour)
inline void RAYLIB_EXPORT convertIntToColour(int x, RGBA &colour)
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
inline int RAYLIB_EXPORT convertColourToInt(const RGBA &colour)
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
}  // namespace ray
#endif  // RAYLIB_RAYEXTRACT_TREES_H