// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFOREST_H
#define RAYLIB_RAYFOREST_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raygrid2d.h"
#include "raylib/rayforeststructure.h"
#include "raylib/raylibconfig.h"
#include "raytreenode.h"
#include "raytrunks.h"

namespace ray
{
/// Class for storing and extracting basic forest information: the location and size of its trees
/// This is the class used in the tool raycloud extract forest, for converting ray clouds to a description
/// of its trees
class RAYLIB_EXPORT Forest
{
public:
  Forest()
    : verbose(true)
    , max_tree_canopy_width(25.0)
    , min_area_(25)
    , undercroft_height(1.5)
    , approx_height_per_radius_(50.0)
    , smooth_iterations_(15)
    , drop_ratio_(0.1)
  {}
  /// Extract tree locations from a specified ray cloud, ground mesh and an optional set of identified trunks
  ForestStructure extract(const std::string &cloud_name_stub, Mesh &mesh,
                          const std::vector<std::pair<Eigen::Vector3d, double>> &trunks, double voxel_width);
  /// Extract tree locations directly from 2D fields: a height field of the canopy (highs) the ground (lows) and the
  /// free space (space).
  ForestStructure extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, const Eigen::ArrayXXd &space,
                          double voxel_width, const std::string &cloud_name_stub);

  /// save the forest data
  bool save(const std::string &filename);

  // parameters
  bool verbose;
  double max_tree_canopy_width;  // maximum width in metres of the widest tree expected. This helps divide up nebulous
                                 // regions
  int min_area_;  // minimum number of pixels that a valid tree should occupy. Below this is considered not enough
                  // evidence and discarded
  double undercroft_height;  // points below this height above the ground are discarded and treated as undercroft
  double approx_height_per_radius_;  // a default ratio for approximate radius estimation from tree height, when trunk
                                     // not observed
  int smooth_iterations_;            // number of iterations for smoothing
  double drop_ratio_;                // a tree edge is one that drops more than this height ratio over one pixel

private:
  /// render the height field to an image file
  void drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield);
  /// find space below the tree node, and output the found space as the @c tip position
  bool findSpace(const TreeNode &node, Eigen::Vector3d &tip);
  /// smooth the height field using repeated averaging of the Moore neighbourhood
  void smoothHeightfield();
  /// include the trunks in the canopy height, as a soft hint for segmentation
  void addTrunkHeights();
  /// render the final segmentation to an image file
  void drawFinalSegmentation(const std::string &filename, std::vector<TreeNode> &trees);
  /// render the hierarchical watershed result
  void renderWatershed(const std::string &cloud_name_stub, std::vector<TreeNode> &trees, std::set<int> &indices);
  /// perform the hierarchical watershed algorithm to segment the trees based on convex crown shapes
  void hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads);

  double voxel_width_;
  Eigen::ArrayXXd heightfield_;
  Eigen::ArrayXXd original_heightfield_;
  Eigen::ArrayXXd lowfield_;
  Eigen::ArrayXXd spacefield_;
  std::vector<std::pair<Eigen::Vector3d, double>> trunks_;
  Eigen::Vector3d min_bounds_, max_bounds_;
  Eigen::ArrayXXi indexfield_;
};

}  // namespace ray

#endif  // RAYLIB_RAYFOREST_H
