// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayconvexhull.h"

#if RAYLIB_WITH_QHULL

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullPointSet.h>
#include <libqhullcpp/QhullPoints.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullVertexSet.h>

#include <map>

namespace ray
{
class lessThan
{
public:
  bool operator()(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) const
  {
    return v1[0] + kPi * v1[1] + kPi * kPi * v1[2] < v2[0] + kPi * v2[1] + kPi * kPi * v2[2];
  }
};

void ConvexHull::construct(const std::vector<Eigen::Vector3d> &points, const Eigen::Vector3d ignoreDirection)
{
  if (points.size() < 3)  // two or fewer points generate an empty mesh
  {
    return;
  }
  std::vector<double> coordinates(points.size() * 3);
  for (int i = 0; i < (int)points.size(); i++)
  {
    coordinates[3 * i + 0] = points[i][0];
    coordinates[3 * i + 1] = points[i][1];
    coordinates[3 * i + 2] = points[i][2];
  }

  orgQhull::Qhull hull;
  hull.setOutputStream(&std::cout);
  hull.runQhull("", 3, int(points.size()), coordinates.data(), "Qt");

  orgQhull::QhullFacetList facets = hull.facetList();
  std::cout << "number of triangles: " << facets.size() << std::endl;
  mesh_.indexList().reserve(facets.size());
  std::cout << "ignore direction: " << ignoreDirection.transpose() << std::endl;
  int count = 0;
  for (const orgQhull::QhullFacet &f : facets)
  {
    double *c = f.hyperplane().coordinates();
    Eigen::Vector3d coords(c[0], c[1], c[2]);
    if (coords.dot(ignoreDirection) <= 0.0)
    {
      count++;
      orgQhull::QhullVertexSet verts = f.vertices();
      int i = 0;
      Eigen::Vector3i index;
      for (const orgQhull::QhullVertex &v : verts) index[i++] = v.point().id();
      Eigen::Vector3d norm = (points[index[1]] - points[index[0]]).cross(points[index[2]] - points[index[0]]);
      if (norm.dot(coords) < 0.0)
        std::swap(index[1], index[2]);
      mesh_.indexList().push_back(index);
    }
  }
  std::cout << "num remaining triangles: " << count << std::endl;
}

ConvexHull::ConvexHull(const std::vector<Eigen::Vector3d> &points)
{
  mesh_.vertices() = points;
}

void ConvexHull::growOutwards(double maxCurvature)
{
  Eigen::Vector3d centre = mean(mesh_.vertices());
  std::vector<Eigen::Vector3d> points = mesh_.vertices();
  for (auto &p : points)
  {
    p -= centre;
    p *= std::pow(p.squaredNorm(), (-1.0 / (1.0 + maxCurvature) - 1.0) * 0.5);
  }

  construct(points, Eigen::Vector3d(0, 0, 0));
}

void ConvexHull::growInwards(double maxCurvature)
{
  Eigen::Vector3d centre = mean(mesh_.vertices());
  std::vector<Eigen::Vector3d> points = mesh_.vertices();
  for (auto &p : points)
  {
    p -= centre;
    p *= std::pow(p.squaredNorm(), (1.0 / (1.0 + maxCurvature) - 1.0) * 0.5);
  }

  construct(points, Eigen::Vector3d(0, 0, 0));
}

void ConvexHull::growInDirection(double maxCurvature, const Eigen::Vector3d &dir)
{
  Eigen::Vector3d centre = mean(mesh_.vertices());
  std::vector<Eigen::Vector3d> points = mesh_.vertices();

  for (auto &p : points)
  {
    Eigen::Vector3d flat = p - centre;
    flat -= dir * dir.dot(flat);
    p += dir * 0.5 * flat.squaredNorm() *
         maxCurvature;  // 0.5 * curv * x^2 means the second differential (the curvature) w.r.t. x is curv
  }

  construct(points, dir);
}
}  // namespace ray
#endif