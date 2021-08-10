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

#define TEST_SCALE 1.0

namespace ray
{
/// For storing and extracting basic forest information 
class RAYLIB_EXPORT Forest
{
public:
  Forest() : verbose(true), undercroft_height(1.0)
  {
  }
  bool extract(const std::string &cloud_name, Mesh &mesh);
  void extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, double voxel_width);
  struct Result
  {
    Eigen::Vector3d tree_tip; // this is in units of pixels horizontally, and metres vertically!
    double ground_height;
    double radius, curvature;
  };

  // in rayforest_draw.cpp
  void drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield);
  void drawTrees(const std::string &filename, const std::vector<Forest::Result> &results, int width, int height);
  bool save(const std::string &filename);

  // parameters
  bool verbose;
  double undercroft_height;
  double tree_roundness;

private:
  double voxel_width_;
  Eigen::ArrayXXd heightfield_;
  Eigen::ArrayXXd lowfield_;
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

}  // namespace ray

#endif  // RAYLIB_RAYFOREST_H
