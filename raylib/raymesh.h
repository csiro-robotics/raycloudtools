// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYMESH_H
#define RAYLIB_RAYMESH_H

#include "raylib/raylibconfig.h"

#include "raycloud.h"
#include "rayutils.h"
#include "rayunused.h"

namespace ray
{
/// A triangular mesh data structure. For mesh based operations.
class RAYLIB_EXPORT Mesh
{
public:
  /// Use the mesh to split a @c cloud based on which side of the mesh its end points are on
  /// The two resulting clouds are @c inside and @c outside
  bool splitCloud(const std::string &cloud_name, double offset, const std::string &inside_name, const std::string &outside_name);

  /// Convert the mesh into a height field (2D array of heights) based on the supplied bounding box and cell width
  void toHeightField(Eigen::ArrayXXd &field, const Eigen::Vector3d &box_min, Eigen::Vector3d box_max,
                     double width, bool fill_gaps = true) const;

  /// access the mesh's vertices
  inline std::vector<Eigen::Vector3d> &vertices() { return vertices_; }
  inline const std::vector<Eigen::Vector3d> &vertices() const { return vertices_; }
  /// access the mesh's index list
  inline std::vector<Eigen::Vector3i> &indexList() { return index_list_; }
  inline const std::vector<Eigen::Vector3i> &indexList() const { return index_list_; }
  /// access the mesh's uv (texture coordinate) list
  inline std::vector<Eigen::Vector3cf> &uvList() { return uv_list_; }
  inline const std::vector<Eigen::Vector3cf> &uvList() const { return uv_list_; }
  /// access the mesh's colours
  inline std::vector<RGBA> &colours() { return colours_; }
  inline const std::vector<RGBA> &colours() const { return colours_; }
  // access the mesh's optional texture name
  inline std::string &textureName() { return texture_name_; }
  inline const std::string &textureName() const { return texture_name_; }

  /// Get first and second order moments of mesh. This can be used as a simple way to compare meshes
  /// numerically. Note that different stats guarantee different meshes, but same stats do not guarantee same meshes
  /// These stats are arranged as the mean vertex location, then the standard deviation in each axis
  Eigen::Array<double, 6, 1> getMoments() const;

  // remove surplus points that are not part of any triangles
  void reduce();

  void translate(const Eigen::Vector3d &offset)
  {
    for (auto &vert: vertices_)
    {
      vert += offset;
    }
  }

private:
  std::vector<Eigen::Vector3d> vertices_;
  std::vector<Eigen::Vector3i> index_list_; // one per triangle, gives the index into the vertices_ array for each corner

  // Optional attributes
  std::vector<Eigen::Vector3cf> uv_list_; // one complex value per face vertex. This is optional
  std::vector<RGBA> colours_;  // optional, if empty then not used
  std::string texture_name_; // optional texture to use if UVs are specifies
};

class Triangle
{
public:
  Eigen::Vector3d corners[3];
  Eigen::Vector3d normal;
  bool tested;
  bool intersectsRay(const Eigen::Vector3d &ray_start, const Eigen::Vector3d &ray_end, double &depth)
  {
    // 1. plane test:
    double d1 = (ray_start - corners[0]).dot(normal);
    double d2 = (ray_end - corners[0]).dot(normal);
    if (d1 * d2 > 0.0)
      return false;

    depth = d1 / (d1 - d2);
    Eigen::Vector3d contact_point = ray_start + (ray_end - ray_start) * depth;

    // next we have to test every sideways direction
    for (int i = 0; i < 3; i++)
    {
      Eigen::Vector3d side = (corners[(i + 1) % 3] - corners[i]).cross(normal);
      if ((contact_point - corners[i]).dot(side) >= 0.0)
        return false;
    }
    return true;
  }
  double distSqrToPoint(const Eigen::Vector3d &point)
  {
    Eigen::Vector3d pos = point - normal * (point - corners[0]).dot(normal);
    bool outs[3];
    double ds[3];
    Eigen::Vector3d sides[3];
    for (int i = 0; i < 3; i++)
    {
      sides[i] = (corners[(i + 1) % 3] - corners[i]).cross(normal);
      ds[i] = (pos - corners[i]).dot(sides[i]);
      outs[i] = ds[i] > 0.0;
    }
    if (outs[0] && outs[1])
      pos = corners[1];
    else if (outs[1] && outs[2])
      pos = corners[2];
    else if (outs[2] && outs[0])
      pos = corners[0];
    else if (outs[0])
      pos -= sides[0] * ds[0] / sides[0].squaredNorm();
    else if (outs[1])
      pos -= sides[1] * ds[1] / sides[1].squaredNorm();
    else if (outs[2])
      pos -= sides[2] * ds[2] / sides[2].squaredNorm();
    return (point - pos).squaredNorm();
  }
  bool intersectsCube(const Eigen::Vector3d &cube_min, double cube_width)
  {
    RAYLIB_UNUSED(cube_min);
    RAYLIB_UNUSED(cube_width);
    // TODO: fill in
    return true;
  }
};
}  // namespace ray

#endif  // RAYLIB_RAYMESH_H
