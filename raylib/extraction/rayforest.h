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
struct RAYLIB_EXPORT TreeNode; 

/// For storing and extracting basic forest information 
class RAYLIB_EXPORT Forest
{
public:
  Forest() : max_tree_canopy_width(25), maximum_drop_within_tree(2.5e10), undercroft_height(1.5)
  {
  }
  void extract(const Cloud &cloud, Mesh &mesh);
  void extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, double voxel_width);
  struct Result
  {
    Eigen::Vector3d tree_tip;
    double ground_height;
    double radius, curvature;
  };

  // in rayforest_draw.cpp
  void drawLowfield(const std::string &filename, const std::vector<TreeNode> &trees);
  void drawSegmentation(const std::string &filename, std::vector<TreeNode> &trees);
  void drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield);
  void drawGraph(const std::string &filename, const std::vector<Vector4d> &data, double x_min, double x_max, double y_max, double strength_max, double a, double b);
  void drawTrees(const std::string &filename, const std::vector<Forest::Result> &results, int width, int height);
  void drawTreeShapes(const std::string &filename, const std::vector<TreeNode> &results, int width, int height);
  void drawFinalSegmentation(const std::string &filename, std::vector<TreeNode> &trees, std::vector<int> &indices);

  // parameters
  bool verbose;
  double tree_roundness;
  double max_tree_canopy_width;  
  double maximum_drop_within_tree;
  double undercroft_height;

private:
  void hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads);
  void calculateTreeParaboloids(std::vector<TreeNode> &trees);
  double estimateRoundnessAndGroundHeight(std::vector<TreeNode> &trees);
  double searchTrees(const std::vector<TreeNode> &trees, int ind, double length_per_radius, std::vector<int> &indices);
  double voxel_width_;
  Eigen::ArrayXXd heightfield_;
  Eigen::ArrayXXd lowfield_;
  Eigen::ArrayXXi indexfield_;
  std::vector<Result> results_;
  Eigen::Vector3d min_bounds_, max_bounds_;
};

// A 2D field structure, like a scalar field or vector field
// Normally I would use an Eigen::ArrayXX, but it only works on scalar types, I'm using below on a Col type.
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
    children[0] = children[1] = -1;
    peak.setZero();
    ground_height = 0;
  }
  TreeNode(int i, int j, double height_, double voxel_width) // TODO: should this be x,y or a distance in metres? probably x,y
  {
    attaches_to = -1;
    min_bound = max_bound = Eigen::Vector2i(i,j);
    double x = (double)i * voxel_width;
    double y = (double)j * voxel_width;
    children[0] = children[1] = -1;
    ground_height = 0;
    peak = Eigen::Vector3d(x, y, height_); // in which case peak should probably be in metres horizontally
  }
  // for calculating paraboloid of best fit:
  struct Node
  {
    Node(){ clear(); }
    void clear() { curv_mat.setZero(); curv_vec.setZero(); abcd.setZero(); }
    inline double avgHeight() const { return curv_vec[3] / area(); }
    inline double height() const { return abcd[3] - (abcd[1]*abcd[1] + abcd[2]*abcd[2])/(4*abcd[0]); }
    inline double heightAt(double x, double y) const { return abcd[0]*(x*x + y*y) + abcd[1]*x + abcd[2]*y + abcd[3]; }
    inline Eigen::Vector3d tip() const { return Eigen::Vector3d(-abcd[1]/(2*abcd[0]), -abcd[2]/(2*abcd[0]), height()); }
    inline double curvature() const { return abcd[0]; }
    inline double crownRadius() const { return 1.0 / -abcd[0]; }
    inline double area() const { return curv_mat(3, 3); }
    inline Eigen::Vector3d mean() const { return Eigen::Vector3d(curv_mat(1,3), curv_mat(2,3), curv_vec[3]) / area(); }

    Eigen::Matrix4d curv_mat;
    Eigen::Vector4d curv_vec;
    Eigen::Vector4d abcd; // the solved paraboloid

    inline void add(double x, double y, double z, double weight) // TODO: this should probably be in SI units
    {
      Eigen::Vector4d vec(x*x + y*y, x, y, 1.0); 
      curv_mat += weight * vec * vec.transpose();
      curv_vec += weight * z*vec;
    }
  };
  Node node;
  Eigen::Vector2i min_bound, max_bound;
  Eigen::Vector3d peak;
  double ground_height;
  double approx_radius;
  int attaches_to;
  int children[2];

//  Eigen::Vector2d centroid() const { return Eigen::Vector2d(curv_mat(1,3) / area(), curv_mat(2,3) / area()); }
 // inline Eigen::Vector3d weightedMean() const { return Eigen::Vector3d(curv_vec[1] / curv_vec[3], curv_vec[2] / curv_vec[3], peak[2]); }
  inline bool validParaboloid(double max_tree_width, double voxel_width) const 
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

#endif  // RAYLIB_RAYFOREST_H
