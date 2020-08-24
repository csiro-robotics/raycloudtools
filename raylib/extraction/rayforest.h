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
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;

namespace ray
{
// A 2D field structure, like a scalar field or vector field
template <class T>
struct Field2D
{
  Field2D() {}
  Field2D(int x, int y){ init(Eigen::Vector2i(x, y)); }
  Field2D(const Eigen::Vector2i &dimensions){ init(dimensions); }
  inline void init(const Eigen::Vector2i &dimensions)
  {
    dims = dimensions;
    data.resize(dims[0] * dims[1]);
  }
  inline T &operator()(const Eigen::Vector2i &ind){ return data[ind[0] + dims[0]*ind[1]]; } 
  inline const T &operator()(const Eigen::Vector2i &ind) const { return data[ind[0] + dims[0]*ind[1]]; }
  inline T &operator()(int x, int y){ return data[x + dims[0]*y]; } 
  inline const T &operator()(int x, int y) const { return data[x + dims[0]*y]; }
  std::vector<T> data;
  Eigen::Vector2i dims;
};

struct RAYLIB_EXPORT TreeNode
{
  TreeNode() : min_bound(1e10,1e10), max_bound(-1e10,-1e10), attaches_to(-1) 
  {
    curv_mat.setZero();
    curv_vec.setZero();
    sum_square_residual = 0.0;
    sum_square_total = 0.0;
    children[0] = children[1] = -1;
    peak.setZero();
  }
  TreeNode(int x, int y, double height_)
  {
    curv_mat.setZero();
    curv_vec.setZero();
    attaches_to = -1;
    min_bound = max_bound = Eigen::Vector2i(x,y);
    addSample(x,y,height_);
    sum_square_residual = 0.0;
    sum_square_total = 0.0;
    children[0] = children[1] = -1;
    peak = Eigen::Vector3d((double)x, (double)y, height_);
  }
  // for calculating paraboloid of best fit:
  Eigen::Matrix4d curv_mat;
  Eigen::Vector4d curv_vec;
  Eigen::Vector4d abcd; // the solved paraboloid
  double sum_square_residual;
  double sum_square_total;
  Eigen::Vector2i min_bound, max_bound;
  Eigen::Vector3d peak;
  int attaches_to;
  int children[2];

  Eigen::Vector2d centroid() const { return Eigen::Vector2d(curv_mat(1,3) / area(), curv_mat(2,3) / area()); }
  inline double area() const { return curv_mat(3,3); }
  inline double avgHeight() const { return curv_vec[3] / area(); }
  inline double height() const { return abcd[3] - (abcd[1]*abcd[1] + abcd[2]*abcd[2])/(4*abcd[0]); }
  inline Eigen::Vector3d tip() const { return Eigen::Vector3d(-abcd[1]/(2*abcd[0]), -abcd[2]/(2*abcd[0]), height()); }
  inline double heightAt(double x, double y) const { return abcd[0]*(x*x + y*y) + abcd[1]*x + abcd[2]*y + abcd[3]; }
  inline double crownRadius() const { return 1.0 / -abcd[0]; }
  inline bool validParaboloid(double max_tree_width) const 
  {
    const double minimum_crown_radius = 0.5;
    const double maximum_crown_radius = max_tree_width; // setting radius to the tree diameter (i.e. twice) as it is an outer bound
    double r = crownRadius();
    if (r<minimum_crown_radius || r > maximum_crown_radius)
      return false;
    Eigen::Vector3d top = tip();
    for (int i = 0; i<2; i++)
      if (top[i] < (double)min_bound[i] || top[i] > (double)max_bound[i])
        return false;
    return true;
  }
  inline void addSample(double x, double y, double z)
  {
    Eigen::Vector4d vec(x*x + y*y, x, y, 1.0);
    curv_mat += vec * vec.transpose();
    curv_vec += z*vec;
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

/// For storing and extracting basic forest information 
class RAYLIB_EXPORT Forest
{
public:
  Forest() : tree_roundness(0), average_height(0), max_tree_canopy_width(22.0) {}
  void extract(const Cloud &cloud);
  void extract(const Eigen::ArrayXXd &heights, double voxel_width);
  struct Result
  {
    std::vector<Eigen::Vector3d> tree_tips;
    double ground_height;
    double treelength_per_crownradius;
  };

  // in rayforest_draw.cpp
  void drawSegmentation(const std::string &filename, const std::vector<TreeNode> &trees);
  void drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield);
  void drawGraph(const std::string &filename, const std::vector<Eigen::Vector3d> &data, double x_max, double y_max, double strength_max);
  void drawTrees(const std::string &filename, const Forest::Result &result);

  double tree_roundness;
  double average_height;
  bool verbose;
  double max_tree_canopy_width; 
private:
  void hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads);
  void calculateTreeParaboloids(std::vector<TreeNode> &trees);
  double estimateRoundnessAndGroundHeight(std::vector<TreeNode> &trees);
  void estimateRoundness(const Mesh &mesh);
  void searchTrees(const std::vector<TreeNode> &trees, int ind, double error, double length_per_radius, double ground_height, std::vector<int> &indices);
  double voxel_width_;
  Eigen::ArrayXXd heightfield_;
  Eigen::ArrayXXi indexfield_;
  Result result_;
};



}  // namespace ray

#endif  // RAYLIB_RAYFOREST_H
