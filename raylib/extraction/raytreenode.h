// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREENODE_H
#define RAYLIB_RAYTREENODE_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;

namespace ray
{
/// The individual nodes representing trees within rayforest's watershed algorithm
struct RAYLIB_EXPORT TreeNode
{
  TreeNode()
    : min_bound(std::numeric_limits<double>::max(), std::numeric_limits<double>::max())
    , max_bound(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest())
    , attaches_to(-1)
    , trunk_id(-1)
  {
    children[0] = children[1] = -1;
    peak.setZero();
    ground_height = height = 0;
    area = 0;
  }
  TreeNode(int i, int j, double height_, double voxel_width,
           int trunkid)  // TODO: should this be x,y or a distance in metres? probably x,y
  {
    attaches_to = -1;
    min_bound = max_bound = Eigen::Vector2i(i, j);
    double x = static_cast<double>(i) * voxel_width;
    double y = static_cast<double>(j) * voxel_width;
    children[0] = children[1] = -1;
    ground_height = height = 0;
    peak = Eigen::Vector3d(x, y, height_);  // in which case peak should probably be in metres horizontally
    trunk_id = trunkid;
    area = 0;
  }

  Eigen::Vector2i min_bound, max_bound;
  Eigen::Vector3d peak;  // in units of metres
  double ground_height;
  double height;
  int attaches_to;
  int children[2];
  int trunk_id;
  int area;

  void updateBound(const Eigen::Vector2i &bmin, const Eigen::Vector2i &bmax)
  {
    for (int i = 0; i < 2; i++)
    {
      min_bound[i] = std::min(min_bound[i], bmin[i]);
      max_bound[i] = std::max(max_bound[i], bmax[i]);
    }
  }
};


}  // namespace ray

#endif  // RAYLIB_RAYTREENODE_H
