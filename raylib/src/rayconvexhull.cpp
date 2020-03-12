// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayconvexhull.h"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullPoints.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullPointSet.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <map>
using namespace std;
using namespace Eigen;
using namespace RAY;

struct lessThan 
{
  bool operator()(const Vector3d &v1, const Vector3d &v2) const { return v1[0] + pi*v1[1] + pi*pi*v1[2] < v2[0] + pi*v2[1] + pi*pi*v2[2]; }
};

void ConvexHull::construct(const vector<Vector3d> &points, const Vector3d ignoreDirection)
{
  vector<double> coordinates(points.size() * 3);
  for (int i = 0; i < (int)points.size(); i++) 
  {
    coordinates[3*i + 0] = points[i][0];
    coordinates[3*i + 1] = points[i][1];
    coordinates[3*i + 2] = points[i][2];
  }

  orgQhull::Qhull hull;
  hull.setOutputStream(&cout);
  hull.runQhull("", 3, points.size(), coordinates.data(), "Qt");

  orgQhull::QhullFacetList facets = hull.facetList();
  cout << "number of triangles: " << facets.size() << endl;
  mesh.indexList.reserve(facets.size());
  cout << "ignore direction: " << ignoreDirection.transpose() << endl;
  int count = 0;
  for (const orgQhull::QhullFacet &f : facets) 
  {
    double *c = f.hyperplane().coordinates();
    Vector3d coords(c[0], c[1], c[2]);
    if (coords.dot(ignoreDirection) <= 0.0)
    {
      count++;
      orgQhull::QhullVertexSet verts = f.vertices();
      int i = 0;
      Vector3i index;
      for (const orgQhull::QhullVertex &v : verts) 
        index[i++] = v.point().id();
      Vector3d norm = (points[index[1]] - points[index[0]]).cross(points[index[2]]-points[index[0]]);
      if (norm.dot(coords) < 0.0)
        swap(index[1], index[2]);
      mesh.indexList.push_back(index);
    }
  }
  cout << "num remaining triangles: " << count << endl;
}

ConvexHull::ConvexHull(const vector<Vector3d> &points)
{
  mesh.vertices = points;
}

void ConvexHull::growOutwards(double maxCurvature)
{
  Vector3d centre = mean(mesh.vertices);
  vector<Vector3d> points = mesh.vertices;
  for (auto &p: points)
  {
    p -= centre;
    p *= pow(p.squaredNorm(), (-1.0/(1.0 + maxCurvature) - 1.0)*0.5);
  }

  construct(points, Vector3d(0,0,0));
}

void ConvexHull::growInwards(double maxCurvature)
{
  Vector3d centre = mean(mesh.vertices);
  vector<Vector3d> points = mesh.vertices;
  for (auto &p: points)
  {
    p -= centre;
    p *= pow(p.squaredNorm(), (1.0/(1.0 + maxCurvature) - 1.0)*0.5);
  }

  construct(points, Vector3d(0,0,0));
}

void ConvexHull::growInDirection(double maxCurvature, const Vector3d &dir)
{
  Vector3d centre = mean(mesh.vertices);
  vector<Vector3d> points = mesh.vertices;

  for (auto &p: points)
  {
    Vector3d flat = p - centre;
    flat -= dir*dir.dot(flat);
    p += dir*flat.squaredNorm()*maxCurvature;
  }

  construct(points, dir);
}

