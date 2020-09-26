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

  /// access the mesh's vertices
  inline std::vector<Eigen::Vector3d> &vertices(){ return vertices_; }
  inline const std::vector<Eigen::Vector3d> &vertices() const { return vertices_; }
  /// access the mesh's index list
  inline std::vector<Eigen::Vector3i> &index_list(){ return index_list_; }
  inline const std::vector<Eigen::Vector3i> &index_list() const { return index_list_; }

  Eigen::Array<double, 6, 1> getMoments() const;
private:
  std::vector<Eigen::Vector3d> vertices_;
  std::vector<Eigen::Vector3i> index_list_; // one per triangle, gives the index into the vertices_ array for each corner
};


}  // namespace ray

#endif  // RAYLIB_RAYMESH_H
