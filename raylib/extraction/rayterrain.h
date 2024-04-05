// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTERRAIN_H
#define RAYLIB_RAYTERRAIN_H

#include "../rayalignment.h"
#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"

typedef Eigen::Matrix<double, 4, 1> Vector4d;

namespace ray
{
/// A class for storage and extraction of a smooth ground height function, as a triangle mesh
/// This is the class used in the tool rayextract terrain, which creates a triangle mesh of the ray cloud's
/// ground surface
class RAYLIB_EXPORT Terrain
{
public:
  /// Extracts a robust smooth surface mesh from the ray cloud. This is a highest lower bound
  /// conditioned on a maximum ground gradient @c gradient. As such it is a sand model of ground extraction,
  /// it treats ground like sand, being incapable of having a gradient beyond the specified value
  /// The input is the @c cloud and its @c file_prefix (to name debug outputs), and a specified @c gradient
  /// The output is the stored mesh, which is accessed with the mesh() accessor.
  void extract(const Cloud &cloud, const Eigen::Vector3d &offset, const std::string &file_prefix, double gradient, bool verbose);

  /// Direct extraction of the pareto front points
  void growUpwards(const std::vector<Eigen::Vector3d> &positions, double gradient);
  void growDownwards(const std::vector<Eigen::Vector3d> &positions, double gradient);

  /// performs voxel-based culling prior to growing upwards
  void growUpwardsFast(const std::vector<Eigen::Vector3d> &ends, double pixel_width, const Eigen::Vector3d &min_bound,
                       const Eigen::Vector3d &max_bound, double gradient);

  /// access the generated mesh
  Mesh &mesh() { return mesh_; }
  const Mesh &mesh() const { return mesh_; }

private:
  Mesh mesh_;
  static void getParetoFront(const std::vector<Vector4d> &points, std::vector<Vector4d> &front);
};


}  // namespace ray

#endif  // RAYLIB_RAYTERRAIN_H
