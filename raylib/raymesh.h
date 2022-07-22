// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYMESH_H
#define RAYLIB_RAYMESH_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"

namespace ray
{
/// A triangular mesh data structure. For mesh based operations. 
class RAYLIB_EXPORT Mesh
{
public:
  /// Use the mesh to split a @c cloud based on which side of the mesh its end points are on
  /// The two resulting clouds are @c inside and @c outside
  void splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside);

  /// Convert the mesh into a height field (2D array of heights) based on the supplied bounding box and cell width
  void toHeightField(Eigen::ArrayXXd &field, const Eigen::Vector3d &box_min, Eigen::Vector3d box_max, double width) const;

  /// access the mesh's vertices
  inline std::vector<Eigen::Vector3d> &vertices(){ return vertices_; }
  inline const std::vector<Eigen::Vector3d> &vertices() const { return vertices_; }
  /// access the mesh's index list
  inline std::vector<Eigen::Vector3i> &indexList(){ return index_list_; }
  inline const std::vector<Eigen::Vector3i> &indexList() const { return index_list_; }
  /// access the mesh's colours
  inline std::vector<RGBA> &colours(){ return colours_; }
  inline const std::vector<RGBA> &colours() const { return colours_; }

  /// Get first and second order moments of mesh. This can be used as a simple way to compare meshes
  /// numerically. Note that different stats guarantee different meshes, but same stats do not guarantee same meshes
  /// These stats are arranged as the mean vertex location, then the standard deviation in each axis
  Eigen::Array<double, 6, 1> getMoments() const;

  // remove surplus points that are not part of any triangles
  void reduce();
private:
  std::vector<Eigen::Vector3d> vertices_;
  std::vector<Eigen::Vector3i> index_list_; // one per triangle, gives the index into the vertices_ array for each corner
  std::vector<RGBA> colours_; // optional, if empty then not used
};


}  // namespace ray

#endif  // RAYLIB_RAYMESH_H
