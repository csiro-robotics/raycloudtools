// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayconcavehull.h"

#include "raydebugdraw.h"

#if RAYLIB_WITH_QHULL
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullPoints.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullPointSet.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <map>
#include <unordered_map>

using namespace std;
using namespace ray;
using namespace Eigen;

#ifdef __unix__
#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>

bool kbhit(void)
{
  struct timeval tv;
  fd_set rdfs;

  tv.tv_sec = 0;
  tv.tv_usec = 0;

  FD_ZERO(&rdfs);
  FD_SET(STDIN_FILENO, &rdfs);

  select(STDIN_FILENO + 1, &rdfs, NULL, NULL, &tv);
  return (bool)FD_ISSET(STDIN_FILENO, &rdfs);
}
#else  // windows?
#include <conio.h>
#endif

static const double deadFace = 1e10;

class Hasher
{
public:
  size_t operator()(const Eigen::Vector2i &key) const { return key[0] + key[1]; }
};

ConcaveHull::ConcaveHull(const std::vector<Eigen::Vector3d> &points)
{
  centre = mean(points);
  unordered_map<Eigen::Vector2i, int, Hasher> edgeLookup(points.size() * 2);

  std::cout << "number of points: " << points.size() << std::endl;
  vertices = points;
  vertex_on_surface.resize(vertices.size());
  for (int i = 0; i < (int)vertex_on_surface.size(); i++) vertex_on_surface[i] = false;

  std::vector<double> coordinates(points.size() * 3);
  for (int i = 0; i < (int)points.size(); i++)
  {
    coordinates[3 * i + 0] = points[i][0];
    coordinates[3 * i + 1] = points[i][1];
    coordinates[3 * i + 2] = points[i][2];
  }

  orgQhull::Qhull hull;
  hull.setOutputStream(&std::cout);
  hull.runQhull("", 3, points.size(), coordinates.data(), "d Qbb Qt");

  orgQhull::QhullFacetList facets = hull.facetList();
  int maxFacets = 0;
  for (const orgQhull::QhullFacet &f : facets) maxFacets = max(maxFacets, f.id() + 1);
  std::cout << "number of total facets: " << facets.size() << std::endl;
  tetrahedra.resize(maxFacets);

  int maxTris = 0;
  for (const orgQhull::QhullFacet &f : facets)
  {
    if (f.isUpperDelaunay())
      continue;
    qh_makeridges(hull.qh(), f.getFacetT());
    for (const orgQhull::QhullRidge &r : f.ridges()) maxTris = max(maxTris, r.id() + 1);
  }
  triangles.resize(maxTris);
  std::cout << "maximum number of triangles: " << maxTris << std::endl;

  int c = 0;
  for (const orgQhull::QhullFacet &f : facets)
  {
    if (f.isUpperDelaunay())
      continue;
    int i = 0;
    Tetrahedron tetra;
    tetra.id = f.id();
    orgQhull::QhullVertexSet verts = f.vertices();
    for (const orgQhull::QhullVertex &v : verts)
    {
      const double *data = v.point().coordinates();
      Eigen::Vector3d vec(data[0], data[1], data[2]);
      tetra.vertices[i++] = v.point().id();
      if ((vertices[v.point().id()] - vec).squaredNorm() > 1e-8)
        std::cout << "vertex data doesn't match its id" << std::endl;
    }

    orgQhull::QhullRidgeSet ridges = f.ridges();
    if (ridges.size() != 4)
      std::cout << "bad number of ridges: " << ridges.size() << std::endl;
    i = 0;
    for (const orgQhull::QhullRidge &r : ridges)
    {
      if (r.vertices().size() != 3)
        std::cout << "bad ridge size: " << r.vertices().size() << std::endl;
      tetra.neighbours[i] = r.topFacet() == f ? r.bottomFacet().id() : r.topFacet().id();
      tetra.triangles[i++] = r.id();
      if (r.id() > maxTris)
        std::cout << "bad rid" << std::endl;
      int j = 0;
      int rid = r.id();
      if (rid >= (int)triangles.size() || r.id() < 0)
        std::cout << "bag bad" << std::endl;
      if (triangles[rid].valid())
      {
        if (triangles[rid].tetrahedra[0] != r.topFacet().id() && triangles[rid].tetrahedra[0] != r.bottomFacet().id())
          std::cout << "replacing with bad data" << std::endl;
      }
      else
      {
        triangles[rid].tetrahedra[0] = r.topFacet().id();
        triangles[rid].tetrahedra[1] = r.bottomFacet().id();
        triangles[rid].is_surface = r.topFacet().isUpperDelaunay() != r.bottomFacet().isUpperDelaunay();
        if (triangles[rid].tetrahedra[0] != f.id() && triangles[rid].tetrahedra[1] != f.id())
          std::cout << "bad too" << std::endl;
        for (const orgQhull::QhullVertex &v : r.vertices())
        {
          int vid = v.point().id();
          if (vid != tetra.vertices[0] && vid != tetra.vertices[1] && vid != tetra.vertices[2] &&
              vid != tetra.vertices[3])
            std::cout << "bad" << std::endl;
          triangles[rid].vertices[j++] = vid;
        }
        for (int i = 0; i < 3; i++)
        {
          int a = triangles[rid].vertices[i];
          int b = triangles[rid].vertices[(i + 1) % 3];
          Eigen::Vector2i v(min(a, b), max(a, b));
          const auto &res = edgeLookup.find(v);
          if (res == edgeLookup.end())
          {
            triangles[rid].edges[i] = edges.size();
            edgeLookup.insert({ v, edges.size() });
            edges.push_back(Edge(v[0], v[1]));
          }
          else
            triangles[rid].edges[i] = res->second;
        }
      }
    }

    if (f.id() >= (int)tetrahedra.size())
      std::cout << "bad" << std::endl;
    tetrahedra[f.id()] = tetra;

    c++;
  }
  std::cout << "number of tetrahedrons: " << c << std::endl;
}

double ConcaveHull::circumcurvature(const ConcaveHull::Tetrahedron &tetra, int triangleID)
{
  Triangle &triangle = triangles[triangleID];
  Eigen::Vector3d vs[4];
  for (int i = 0; i < 4; i++) vs[i] = vertices[tetra.vertices[i]];

  Eigen::Vector3d m1 = (vs[0] + vs[1]) * 0.5;
  Eigen::Vector3d m2 = (vs[0] + vs[2]) * 0.5;
  Eigen::Vector3d normal = (vs[2] - vs[0]).cross(vs[1] - vs[0]);
  Eigen::Vector3d up = (vs[2] - vs[0]).cross(normal);

  // defines the line midway between vs[0], vs[1] and vs[2]
  double t1 = (m1 - m2).dot(vs[1] - vs[0]) / up.dot(vs[1] - vs[0]);
  Eigen::Vector3d midbase = m2 + up * t1;

  // defines the plane midway between vs[0] and vs[3]
  Eigen::Vector3d m = (vs[3] + vs[0]) * 0.5;
  Eigen::Vector3d dir = vs[3] - vs[0];

  // intersection of line and plane is midway between all four points
  double t = (m - midbase).dot(dir) / normal.dot(dir);
  Eigen::Vector3d circumcentre = midbase + normal * t;
  Eigen::Vector3d centre = (vs[0] + vs[1] + vs[2] + vs[3]) * 0.25;

  // return the distance from a vertex to this circumcentre
  double circumradius = (circumcentre - vs[0]).norm();

  Eigen::Vector3d cs[3];
  for (int i = 0; i < 3; i++) cs[i] = vertices[triangle.vertices[i]];

  Eigen::Vector3d triNormal = (cs[2] - cs[0]).cross(cs[1] - cs[0]);
  double circumcentreSide = (circumcentre - cs[0]).dot(triNormal);
  if (circumcentreSide * (centre - cs[0]).dot(triNormal) > 0.0)
  {
    // we are ready to make this a dead face, but before we do, it is possible that the tetrahedron in question includes
    // an existing surface face in which case, the circumcentre might not in fact be above the current surface.

    int faceIntersects = -1;
    int numFaceIntersects = 0;
    for (int j = 0; j < 4; j++)
    {
      if (tetra.triangles[j] == triangleID)
        continue;
      if (triangles[tetra.triangles[j]].surface_face_cached.triangle != -1)
      {
        numFaceIntersects++;
        faceIntersects = j;
      }
    }
    if (numFaceIntersects == 1)
    {
      int otherFace = tetra.triangles[faceIntersects];
      Triangle tri = triangles[otherFace];
      Eigen::Vector3d ds[3];
      for (int i = 0; i < 3; i++) ds[i] = vertices[tri.vertices[i]];

      Eigen::Vector3d triNormal2 = (ds[2] - ds[0]).cross(ds[1] - ds[0]);
      double circumcentreSide2 = (circumcentre - ds[0]).dot(triNormal2);
      bool aboveOtherFace = circumcentreSide2 * (centre - ds[0]).dot(triNormal2) > 0.0;
      if (aboveOtherFace)
        return deadFace;
    }
    else
      return deadFace;
  }

  return 1.0 / circumradius;
}

static int newTriCount = 0;

bool ConcaveHull::growFront(double maxCurvature)
{
  SurfaceFace face = *surface.begin();
  if (face.curvature == deadFace || face.curvature > maxCurvature)
    return false;
  int vertexI = 0;
  for (int i = 0; i < 4; i++)
  {
    int v = tetrahedra[face.tetrahedron].vertices[i];
    Triangle &tri = triangles[face.triangle];
    if (v != tri.vertices[0] && v != tri.vertices[1] && v != tri.vertices[2])
    {
      vertexI = i;
      break;
    }
  }
  Tetrahedron &tetra = tetrahedra[face.tetrahedron];
  int newVertex = tetra.vertices[vertexI];
  bool intersects = vertex_on_surface[newVertex];
  int faceIntersects = -1;
  int numFaceIntersects = 0;
  const SurfaceFace *faceIntersectTri = NULL;
  for (int j = 0; j < 4; j++)
  {
    if (tetra.triangles[j] == face.triangle)
      continue;
    if (triangles[tetra.triangles[j]].surface_face_cached.triangle != -1)
    {
      numFaceIntersects++;
      faceIntersects = j;
      faceIntersectTri = &triangles[tetra.triangles[j]].surface_face_cached;
    }
  }

  surface.erase(surface.begin());
  if (numFaceIntersects == 1)
  {
    intersects = false;
    int otherVertex = 0;
    for (int i = 0; i < 3; i++)
    {
      int v = triangles[face.triangle].vertices[i];
      Triangle &tri = triangles[tetra.triangles[faceIntersects]];
      if (v != tri.vertices[0] && v != tri.vertices[1] && v != tri.vertices[2])
        otherVertex = v;
    }
    int tri2 = tetra.triangles[(faceIntersects + 1) % 3];
    if (tri2 == face.triangle)
      tri2 = tetra.triangles[(faceIntersects + 2) % 3];
    int v0 = min(otherVertex, newVertex);
    int v1 = max(otherVertex, newVertex);
    for (int i = 0; i < 3; i++)
    {
      Edge &edge = edges[triangles[tri2].edges[i]];
      if (edge.vertices[0] == v0 && edge.vertices[1] == v1 && edge.has_had_face)
        intersects = true;
    }
    if (!intersects)
      surface.erase(*faceIntersectTri);
  }
  if (intersects)
  {
    face.curvature = deadFace;
    triangles[face.triangle].surface_face_cached = face;
    surface.insert(face);  // put it at the back of the queue
    return true;
  }

  newTriCount++;
  for (int i = 0; i < 4; i++)
  {
    if (tetra.triangles[i] == face.triangle || (numFaceIntersects == 1 && i == faceIntersects))
      continue;
    SurfaceFace newFace;

    newFace.tetrahedron = tetra.neighbours[i];
    newFace.triangle = tetra.triangles[i];
    for (int j = 0; j < 3; j++) vertex_on_surface[triangles[newFace.triangle].vertices[j]] = true;
    double grad;
    if (triangles[newFace.triangle].is_surface)
      newFace.curvature = grad = deadFace;
    else
      newFace.curvature = circumcurvature(tetrahedra[newFace.tetrahedron], newFace.triangle);
    triangles[newFace.triangle].surface_face_cached = newFace;
    for (int j = 0; j < 3; j++) edges[triangles[newFace.triangle].edges[j]].has_had_face = true;
    surface.insert(newFace);
  }
  return true;
}

void ConcaveHull::growSurface(double maxCurvature)
{
  do
  {
    if (!(newTriCount % 1600))
    {
      std::cout << "max curvature of structure: " << surface.begin()->curvature << std::endl;
      std::vector<std::vector<Eigen::Vector3d>> tris;
      for (auto &tri : surface)
      {
        std::vector<Eigen::Vector3d> corners(3);
        for (int i = 0; i < 3; i++) corners[i] = vertices[triangles[tri.triangle].vertices[i]];
        tris.push_back(corners);
      }
      // if (DebugDraw *debug = DebugDraw::instance())
      // {
      //   debug->drawTriangles(tris, 0.5);
      // }
      if (kbhit())
        return;
    }
  } while (growFront(maxCurvature));
}

// starting with given tetrahedron, grow it outwards to achieve a maximum curvature
void ConcaveHull::growOutwards(const ConcaveHull::Tetrahedron &tetra, double maxCurvature)
{
  surface.clear();
  for (int i = 0; i < 4; i++)
  {
    SurfaceFace face;
    face.tetrahedron = tetra.neighbours[i];
    face.triangle = tetra.triangles[i];
    if (triangles[face.triangle].is_surface)
      face.curvature = deadFace;
    else
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
    triangles[face.triangle].surface_face_cached = face;
    surface.insert(face);
  }
  growSurface(maxCurvature);
}

void ConcaveHull::growOutwards(double maxCurvature)
{
  bool found = false;
  for (auto &tetra : tetrahedra)
  {
    if (insideTetrahedron(centre, tetra))
    {
      growOutwards(tetra, maxCurvature);
      found = true;
    }
  }
  if (!found)
    std::cout << "could not find a tetrahedron at the point cloud mean location" << std::endl;
}

// starting with the outer (convex) surface mesh, grow inwards up to the maxCurvature value
void ConcaveHull::growInwards(double maxCurvature)
{
  surface.clear();
  // find the surface triangles...
  for (int i = 0; i < (int)triangles.size(); i++)
  {
    Triangle &tri = triangles[i];
    if (tri.is_surface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      face.triangle = i;
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
      triangles[i].surface_face_cached = face;
      surface.insert(face);
    }
  }
  growSurface(maxCurvature);
}

void ConcaveHull::growInDirection(double maxCurvature, const Eigen::Vector3d &dir)
{
  surface.clear();
  // find the surface triangles...
  for (int i = 0; i < (int)triangles.size(); i++)
  {
    Triangle &tri = triangles[i];
    if (tri.is_surface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      if (face.tetrahedron < 0)
        std::cout << "bad face tetrahedron" << std::endl;
      Eigen::Vector3d mid(0, 0, 0);
      for (int j = 0; j < 4; j++) mid += vertices[tetrahedra[face.tetrahedron].vertices[j]] / 4.0;
      Eigen::Vector3d normal = (vertices[tri.vertices[2]] - vertices[tri.vertices[0]])
                          .cross(vertices[tri.vertices[1]] - vertices[tri.vertices[0]]);
      if ((mid - vertices[tri.vertices[0]]).dot(normal) < 0.0)
        normal = -normal;
      if (normal.dot(dir) < 0.0)  // downwards facing
        continue;
      face.triangle = i;
      for (int j = 0; j < 3; j++)
      {
        vertex_on_surface[tri.vertices[j]] = true;
        edges[tri.edges[j]].has_had_face = true;
      }
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
      triangles[i].surface_face_cached = face;
      surface.insert(face);
    }
  }
  growSurface(maxCurvature);
}
#endif
