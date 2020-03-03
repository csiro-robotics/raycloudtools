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
  triangles.reserve(facets.size());
  for (const orgQhull::QhullFacet &f : facets) 
  {
    double *c = f.hyperplane().coordinates();
    Vector3d coords(c[0], c[1], c[2]);
    if (coords.dot(ignoreDirection) <= 0)
    {
      orgQhull::QhullVertexSet verts = f.vertices();
      Triangle triangle;
      int i = 0;
      for (const orgQhull::QhullVertex &v : verts) 
      {
        const double *data = v.point().coordinates();
        triangle.vertices[i++] = Vector3d(data[0], data[1], data[2]);
      }
      triangles.push_back(triangle);
    }
  }
}

ConvexHull::ConvexHull(const vector<Vector3d> &points)
{
  this->points = points;
}

void ConvexHull::growOutwards(double maxCurvature)
{
  Vector3d centre = mean(points);
  for (auto &p: points)
  {
    p -= centre;
    p /= pow((p-centre).squaredNorm(), maxCurvature*0.5);
  }

  construct(points, Vector3d(0,0,0));

  for (auto &tri: triangles)
  {
    for (auto &p: tri.vertices)
    {
      p *= pow(p.squaredNorm(), maxCurvature*0.5);
      p += centre;
    }
  }
}

void ConvexHull::growInwards(double maxCurvature)
{
  Vector3d centre = mean(points);
  for (auto &p: points)
  {
    p -= centre;
    p *= pow(p.squaredNorm(), maxCurvature*0.5);
  }

  construct(points, Vector3d(0,0,0));

  for (auto &tri: triangles)
    for (auto &p: tri.vertices)
    {
      p /= pow(p.squaredNorm(), maxCurvature*0.5);
      p += centre;
    }
}

void ConvexHull::growInDirection(double maxCurvature, const Vector3d &dir)
{
  Vector3d centre = mean(points);

  for (auto &p: points)
  {
    Vector3d flat = p - centre;
    flat -= dir*dir.dot(flat);
    p += dir*flat.squaredNorm()*maxCurvature;
  }

  construct(points, dir);

  for (auto &tri: triangles)
  {
    for (auto &p: tri.vertices)
    {
      Vector3d flat = p - centre;
      flat -= dir*dir.dot(flat)*maxCurvature;
      p -= dir*flat.squaredNorm();
    }
  }
}

