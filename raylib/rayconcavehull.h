// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCONCAVEHULL_H
#define RAYLIB_RAYCONCAVEHULL_H

#include <set>
#include "raylib/raylibconfig.h"
#include "raylib/raymesh.h"
#include "rayutils.h"

#if RAYLIB_WITH_QHULL
namespace ray
{
/// Class for calculating concave hulls. This is a single, simply connected polyhedron with a maximum allowed
/// concave curvature @c maxCurvature. It is a generalisation of a convex hull (@c maxCurvature=0) to any
/// concave curvature. It is used like a vacuum wrapping of a ray cloud, particularly to extract ground terrain.
/// It contains multiple methods for generating the hull, depending on which direction it should grow
class RAYLIB_EXPORT ConcaveHull
{
public:
  /// construct the hull from the ray cloud's end points
  ConcaveHull(const std::vector<Eigen::Vector3d> &points);

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
  class SurfaceFace
  {
  public:
    SurfaceFace() { triangle = -1; }
    int tetrahedron;
    int triangle;
    double curvature;
  };
  class Edge
  {
  public:
    Edge() {}
    Edge(int v1, int v2)
    {
      vertices[0] = v1;
      vertices[1] = v2;
      has_had_face = false;
    }
    int vertices[2];
    bool has_had_face;
  };
  class Triangle
  {
  public:
    Triangle()
    {
      vertices[0] = vertices[1] = vertices[2] = -1;
      edges = Eigen::Vector3i(-1, -1, -1);
      is_surface = false;
      used = false;
    }
    bool valid() { return vertices[0] != -1; }
    bool is_surface;
    bool used;
    Eigen::Vector3i vertices;
    Eigen::Vector3i edges;
    int tetrahedra[2];
    SurfaceFace surface_face_cached;
  };
  class Tetrahedron
  {
  public:
    Tetrahedron()
    {
      vertices[0] = vertices[1] = vertices[2] = vertices[3] = -1;
      seen = false;
      id = -1;
    }
    bool valid() { return vertices[0] != -1; }
    int vertices[4];
    int triangles[4];
    int neighbours[4];
    int id;
    bool seen;
  };
  class FaceComp
  {
  public:
    bool operator()(const SurfaceFace &lhs, const SurfaceFace &rhs) const
    {
      if (lhs.curvature == rhs.curvature)
      {
        if (lhs.triangle == rhs.triangle)
          return lhs.tetrahedron < rhs.tetrahedron;
        return lhs.triangle < rhs.triangle;
      }
      return lhs.curvature < rhs.curvature;
    }
  };

  inline bool insideTetrahedron(const Eigen::Vector3d &pos, const Tetrahedron &tetra)
  {
    Eigen::Vector3d mid(0, 0, 0);
    if (tetra.vertices[0] == -1 || tetra.vertices[1] == -1 || tetra.vertices[2] == -1 ||
        tetra.vertices[3] == -1)  // an outer tetrahedron
      return false;
    for (int j = 0; j < 4; j++) mid += vertices_[tetra.vertices[j]] / 4.0;
    for (int i = 0; i < 4; i++)
    {
      Eigen::Vector3d vs[3];
      for (int j = 0; j < 3; j++) vs[j] = vertices_[triangles_[tetra.triangles[i]].vertices[j]];
      Eigen::Vector3d normal = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
      if ((pos - vs[0]).dot(normal) * (mid - vs[0]).dot(normal) < 0)
        return false;
    }
    return true;
  }
  void growSurface(double maxCurvature);
  bool growFront(double maxCurvature);
  double circumcurvature(const ConcaveHull::Tetrahedron &tetra, int triangleID);
  void growOutwards(const ConcaveHull::Tetrahedron &tetra, double maxCurvature);
  void convertToMesh();

  std::vector<bool> vertex_on_surface_;
  std::vector<Eigen::Vector3d> vertices_;
  std::vector<Edge> edges_;
  std::vector<Triangle> triangles_;
  std::vector<Tetrahedron> tetrahedra_;
  Eigen::Vector3d centre_;
  std::set<SurfaceFace, FaceComp> surface_;
  Mesh mesh_;
};
}  // namespace ray

#endif

#endif  // RAYLIB_RAYCONCAVEHULL_H
