// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFOREST_H
#define RAYLIB_RAYFOREST_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include "../raymesh.h"
#include "raybranches.h"
#include "raywatershed.h"
#include "rayoccupancy2d.h"
#include "raytreesummary.h"

namespace ray
{
struct Cluster
{
  Eigen::Vector3d min_bound, max_bound;
  std::vector<int> ids;
  bool active;
  int trunk_id;
};

/// For storing and extracting basic forest information 
class RAYLIB_EXPORT Forest
{
public:
  Forest() : verbose(true), max_tree_canopy_width(25.0), min_area_(25), undercroft_height(1.5) {}
  std::vector<struct TreeSummary> extract(const std::string &cloud_name_stub, Mesh &mesh, const std::vector<std::pair<Eigen::Vector3d, double> > &trunks);
  std::vector<struct TreeSummary> extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, const Eigen::ArrayXXd &space, double voxel_width, const std::string &cloud_name_stub);

  bool save(const std::string &filename);

  // parameters
  bool verbose;
  double max_tree_canopy_width; 
  int min_area_;
  double tree_roundness;
  double undercroft_height;

private:
  void drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield);
  bool findSpace(const Cluster &cluster, const std::vector<Eigen::Vector3d> &points, Eigen::Vector3d &tip);
  bool findSpace2(const TreeNode &node, Eigen::Vector3d &tip);
  void smoothHeightfield();
  void addTrunkHeights();

  double voxel_width_;
  Eigen::ArrayXXd heightfield_;
  Eigen::ArrayXXd original_heightfield_;
  Eigen::ArrayXXd lowfield_;
  Eigen::ArrayXXd spacefield_;
  std::vector<std::pair<Eigen::Vector3d, double> > trunks_;
  Eigen::Vector3d min_bounds_, max_bounds_;

  // agglomerate:
  void agglomerate(const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> index_list, double min_diameter_per_height, double max_diameter_per_height, std::vector<Cluster> &point_clusters);
  void renderAgglomeration(const std::vector<Cluster> &point_clusters, const std::vector<Eigen::Vector3d> &verts, const std::string &cloud_name_stub);


  // watershed:
 // void drawTrees(const std::string &filename, const std::vector<TreeSummary> &results, int width, int height);
  void drawFinalSegmentation(const std::string &filename, std::vector<TreeNode> &trees);
  void renderWatershed(const std::string &cloud_name_stub, std::vector<TreeNode> &trees, std::set<int> &indices);
  void hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads);
  Eigen::ArrayXXi indexfield_;
};

}  // namespace ray

#endif  // RAYLIB_RAYFOREST_H
