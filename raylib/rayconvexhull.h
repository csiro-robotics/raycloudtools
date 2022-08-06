// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCONVEXHULL_H
#define RAYLIB_RAYCONVEXHULL_H

#include "raylib/raylibconfig.h"

#include <set>
#include "raymesh.h"
#include "rayutils.h"
#if RAYLIB_WITH_QHULL

namespace ray
{
/// This class is poorly named (I must rename it!)
/// It wraps a non-convex mesh onto a specified ray cloud by using a convec hull algorithm on a non-linear
/// transformation of the cloud. The result is a non-convex 'vacuum wrapping' of the ray cloud that is computationally
/// much faster than convex hull, but which does not support overhangs. Multiple wrap directions allow it to be used in
/// a variety of situations.
class RAYLIB_EXPORT ConvexHull
{
public:
  /// construct the hull from the end points of a ray cloud
  ConvexHull(const std::vector<Eigen::Vector3d> &points);

  /// inwards growth is for wrapping an object from the outside, such as a plane
  void growInwards(double maxCurvature);
  /// outwards growth is for wrapping a scene from within it, such as for getting a mesh of a room or a cave
  void growOutwards(double maxCurvature);
  /// growth in a single direction
  void growInDirection(double maxCurvature, const Eigen::Vector3d &dir);
  /// upwards growth is for extracting the ground mesh underneath a (potentially cluttered) scene
  void growUpwards(double maxCurvature) { growInDirection(maxCurvature, Eigen::Vector3d(0, 0, 1)); }
  /// downwards growth is for extracting a 'blanket mesh' of a forest canopy or of building roofs
  void growDownwards(double maxCurvature) { growInDirection(maxCurvature, Eigen::Vector3d(0, 0, -1)); }

  /// access the generated mesh
  Mesh &mesh() { return mesh_; }
  const Mesh &mesh() const { return mesh_; }

private:
  Mesh mesh_;
  void construct(const std::vector<Eigen::Vector3d> &points, const Eigen::Vector3d ignoreDirection);
};
}  // namespace ray
#endif

#endif  // RAYLIB_RAYCONVEXHULL_H
