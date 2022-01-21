// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYWATERSHED_H
#define RAYLIB_RAYWATERSHED_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include "../raymesh.h"
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;

namespace ray
{
struct RAYLIB_EXPORT TreeNode
{
  TreeNode() : min_bound(1e10,1e10), max_bound(-1e10,-1e10), attaches_to(-1), trunk_id(-1) 
  {
    children[0] = children[1] = -1;
    peak.setZero();
    ground_height = height = 0;
    area = 0;
  }
  TreeNode(int i, int j, double height_, double voxel_width, int trunkid) // TODO: should this be x,y or a distance in metres? probably x,y
  {
    attaches_to = -1;
    min_bound = max_bound = Eigen::Vector2i(i,j);
    double x = (double)i * voxel_width;
    double y = (double)j * voxel_width;
    children[0] = children[1] = -1;
    ground_height = height = 0;
    peak = Eigen::Vector3d(x, y, height_); // in which case peak should probably be in metres horizontally
    trunk_id = trunkid;
    area = 0;
  }

  Eigen::Vector2i min_bound, max_bound;
  Eigen::Vector3d peak; // in units of metres
  double ground_height;
  double height;
  int attaches_to;
  int children[2];
  int trunk_id;
  int area;

//  Eigen::Vector2d centroid() const { return Eigen::Vector2d(curv_mat(1,3) / area(), curv_mat(2,3) / area()); }
 // inline Eigen::Vector3d weightedMean() const { return Eigen::Vector3d(curv_vec[1] / curv_vec[3], curv_vec[2] / curv_vec[3], peak[2]); }
  inline bool validParaboloid(double /*max_tree_width*/, double /*voxel_width*/) const 
  {
    // Add voxel_width*voxel_width* to below two lines, to verify voxel_width independence
/*    const double minimum_crown_radius = 0.5;
    const double maximum_crown_radius = max_tree_width; // setting radius to the tree diameter (i.e. twice) as it is an outer bound
    double r = node.crownRadius();
    if (r<minimum_crown_radius || r > maximum_crown_radius)
      return false;
    Eigen::Vector3d top = node.tip();
    for (int i = 0; i<2; i++)
      if (top[i] < (double)min_bound[i]*voxel_width || top[i] > (double)max_bound[i]*voxel_width)
        return false;  */
    return true;
  }
  void updateBound(const Eigen::Vector2i &bmin, const Eigen::Vector2i &bmax)
  {
    for (int i = 0; i<2; i++)
    {
      min_bound[i] = std::min(min_bound[i], bmin[i]);
      max_bound[i] = std::max(max_bound[i], bmax[i]);
    }
  }
};



}  // namespace ray

#endif  // RAYLIB_RAYWATERSHED_H
