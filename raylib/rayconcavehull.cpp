// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayconcavehull.h"

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
#include <unordered_map>

#ifdef __unix__
#include <cstdio>
#include <sys/time.h>
#include <sys/types.h>
#include <termios.h>
#include <unistd.h>
#else  // windows?
#include <conio.h>
#endif

namespace ray
{
static const double deadFace = 1e10;

class Hasher
{
public:
  size_t operator()(const Eigen::Vector2i &key) const { return key[0] + key[1]; }
};

ConcaveHull::ConcaveHull(const std::vector<Eigen::Vector3d> &points)
{
  centre_ = mean(points);
  std::unordered_map<Eigen::Vector2i, int, Hasher> edgeLookup(points.size() * 2);

  std::cout << "number of points: " << points.size() << std::endl;
  vertices_ = points;
  vertex_on_surface_.resize(vertices_.size());
  for (int i = 0; i < (int)vertex_on_surface_.size(); i++) vertex_on_surface_[i] = false;

  std::vector<double> coordinates(points.size() * 3);
  for (int i = 0; i < (int)points.size(); i++)
  {
    coordinates[3 * i + 0] = points[i][0];
    coordinates[3 * i + 1] = points[i][1];
    coordinates[3 * i + 2] = points[i][2];
  }

  orgQhull::Qhull hull;
  hull.setOutputStream(&std::cout);
  hull.runQhull("", 3, int(points.size()), coordinates.data(), "d Qbb Qt");

  orgQhull::QhullFacetList facets = hull.facetList();
  int maxFacets = 0;
  for (const orgQhull::QhullFacet &f : facets) maxFacets = std::max(maxFacets, f.id() + 1);
  std::cout << "number of total facets: " << facets.size() << std::endl;
  tetrahedra_.resize(maxFacets);

  int maxTris = 0;
  for (const orgQhull::QhullFacet &f : facets)
  {
    if (f.isUpperDelaunay())
      continue;
    qh_makeridges(hull.qh(), f.getFacetT());
    for (const orgQhull::QhullRidge &r : f.ridges()) maxTris = std::max(maxTris, r.id() + 1);
  }
  triangles_.resize(maxTris);
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
      if ((vertices_[v.point().id()] - vec).squaredNorm() > 1e-8)
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
      if (rid >= (int)triangles_.size() || r.id() < 0)
        std::cout << "bag bad" << std::endl;
      if (triangles_[rid].valid())
      {
        if (triangles_[rid].tetrahedra[0] != r.topFacet().id() && triangles_[rid].tetrahedra[0] != r.bottomFacet().id())
          std::cout << "replacing with bad data" << std::endl;
      }
      else
      {
        triangles_[rid].tetrahedra[0] = r.topFacet().id();
        triangles_[rid].tetrahedra[1] = r.bottomFacet().id();
        triangles_[rid].is_surface = r.topFacet().isUpperDelaunay() != r.bottomFacet().isUpperDelaunay();
        if (triangles_[rid].tetrahedra[0] != f.id() && triangles_[rid].tetrahedra[1] != f.id())
          std::cout << "bad too" << std::endl;
        for (const orgQhull::QhullVertex &v : r.vertices())
        {
          int vid = v.point().id();
          if (vid != tetra.vertices[0] && vid != tetra.vertices[1] && vid != tetra.vertices[2] &&
              vid != tetra.vertices[3])
            std::cout << "bad" << std::endl;
          triangles_[rid].vertices[j++] = vid;
        }
        for (int i = 0; i < 3; i++)
        {
          int a = triangles_[rid].vertices[i];
          int b = triangles_[rid].vertices[(i + 1) % 3];
          Eigen::Vector2i v(std::min(a, b), std::max(a, b));
          const auto &res = edgeLookup.find(v);
          if (res == edgeLookup.end())
          {
            triangles_[rid].edges[i] = int(edges_.size());
            edgeLookup.insert({ v, edges_.size() });
            edges_.push_back(Edge(v[0], v[1]));
          }
          else
            triangles_[rid].edges[i] = res->second;
        }
      }
    }

    if (f.id() >= (int)tetrahedra_.size())
      std::cout << "bad" << std::endl;
    tetrahedra_[f.id()] = tetra;

    c++;
  }
  std::cout << "number of tetrahedrons: " << c << std::endl;
}

double ConcaveHull::circumcurvature(const ConcaveHull::Tetrahedron &tetra, int triangleID)
{
  Triangle &triangle = triangles_[triangleID];
  Eigen::Vector3d vs[4];
  for (int i = 0; i < 4; i++) vs[i] = vertices_[tetra.vertices[i]];

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
  Eigen::Vector3d circumcentre_ = midbase + normal * t;
  Eigen::Vector3d centre_ = (vs[0] + vs[1] + vs[2] + vs[3]) * 0.25;

  // return the distance from a vertex to this circumcentre_
  double circumradius = (circumcentre_ - vs[0]).norm();

  Eigen::Vector3d cs[3];
  for (int i = 0; i < 3; i++) cs[i] = vertices_[triangle.vertices[i]];

  Eigen::Vector3d triNormal = (cs[2] - cs[0]).cross(cs[1] - cs[0]);
  double circumcentre_Side = (circumcentre_ - cs[0]).dot(triNormal);
  if (circumcentre_Side * (centre_ - cs[0]).dot(triNormal) > 0.0)
  {
    // we are ready to make this a dead face, but before we do, it is possible that the tetrahedron in question includes
    // an existing surface face in which case, the circumcentre_ might not in fact be above the current surface.

    int faceIntersects = -1;
    int numFaceIntersects = 0;
    for (int j = 0; j < 4; j++)
    {
      if (tetra.triangles[j] == triangleID)
        continue;
      if (triangles_[tetra.triangles[j]].surface_face_cached.triangle != -1)
      {
        numFaceIntersects++;
        faceIntersects = j;
      }
    }
    if (numFaceIntersects == 1)
    {
      int otherFace = tetra.triangles[faceIntersects];
      Triangle tri = triangles_[otherFace];
      Eigen::Vector3d ds[3];
      for (int i = 0; i < 3; i++) ds[i] = vertices_[tri.vertices[i]];

      Eigen::Vector3d triNormal2 = (ds[2] - ds[0]).cross(ds[1] - ds[0]);
      double circumcentre_Side2 = (circumcentre_ - ds[0]).dot(triNormal2);
      bool aboveOtherFace = circumcentre_Side2 * (centre_ - ds[0]).dot(triNormal2) > 0.0;
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
  SurfaceFace face = *surface_.begin();
  if (face.curvature == deadFace || face.curvature > maxCurvature)
    return false;
  int vertexI = 0;
  for (int i = 0; i < 4; i++)
  {
    int v = tetrahedra_[face.tetrahedron].vertices[i];
    Triangle &tri = triangles_[face.triangle];
    if (v != tri.vertices[0] && v != tri.vertices[1] && v != tri.vertices[2])
    {
      vertexI = i;
      break;
    }
  }
  Tetrahedron &tetra = tetrahedra_[face.tetrahedron];
  int newVertex = tetra.vertices[vertexI];
  bool intersects = vertex_on_surface_[newVertex];
  int faceIntersects = -1;
  int numFaceIntersects = 0;
  const SurfaceFace *faceIntersectTri = nullptr;
  for (int j = 0; j < 4; j++)
  {
    if (tetra.triangles[j] == face.triangle)
      continue;
    if (triangles_[tetra.triangles[j]].surface_face_cached.triangle != -1)
    {
      numFaceIntersects++;
      faceIntersects = j;
      faceIntersectTri = &triangles_[tetra.triangles[j]].surface_face_cached;
    }
  }

  surface_.erase(surface_.begin());
  if (numFaceIntersects == 1)
  {
    intersects = false;
    int otherVertex = 0;
    for (int i = 0; i < 3; i++)
    {
      int v = triangles_[face.triangle].vertices[i];
      Triangle &tri = triangles_[tetra.triangles[faceIntersects]];
      if (v != tri.vertices[0] && v != tri.vertices[1] && v != tri.vertices[2])
        otherVertex = v;
    }
    int tri2 = tetra.triangles[(faceIntersects + 1) % 3];
    if (tri2 == face.triangle)
      tri2 = tetra.triangles[(faceIntersects + 2) % 3];
    int v0 = std::min(otherVertex, newVertex);
    int v1 = std::max(otherVertex, newVertex);
    for (int i = 0; i < 3; i++)
    {
      Edge &edge = edges_[triangles_[tri2].edges[i]];
      if (edge.vertices[0] == v0 && edge.vertices[1] == v1 && edge.has_had_face)
        intersects = true;
    }
    if (!intersects)
      surface_.erase(*faceIntersectTri);
  }
  if (intersects)
  {
    face.curvature = deadFace;
    triangles_[face.triangle].surface_face_cached = face;
    surface_.insert(face);  // put it at the back of the queue
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
    for (int j = 0; j < 3; j++) vertex_on_surface_[triangles_[newFace.triangle].vertices[j]] = true;
    double grad;
    if (triangles_[newFace.triangle].is_surface)
      newFace.curvature = grad = deadFace;
    else
      newFace.curvature = circumcurvature(tetrahedra_[newFace.tetrahedron], newFace.triangle);
    triangles_[newFace.triangle].surface_face_cached = newFace;
    for (int j = 0; j < 3; j++) edges_[triangles_[newFace.triangle].edges[j]].has_had_face = true;
    surface_.insert(newFace);
  }
  return true;
}

void ConcaveHull::growSurface(double maxCurvature)
{
  do
  {
    if (!(newTriCount % 1600))
    {
      std::cout << "max curvature of structure: " << surface_.begin()->curvature << std::endl;
      std::vector<std::vector<Eigen::Vector3d>> tris;
      for (auto &tri : surface_)
      {
        std::vector<Eigen::Vector3d> corners(3);
        for (int i = 0; i < 3; i++) corners[i] = vertices_[triangles_[tri.triangle].vertices[i]];
        tris.push_back(corners);
      }
    }
  } while (growFront(maxCurvature));
}

// starting with given tetrahedron, grow it outwards to achieve a maximum curvature
void ConcaveHull::growOutwards(const ConcaveHull::Tetrahedron &tetra, double maxCurvature)
{
  surface_.clear();
  for (int i = 0; i < 4; i++)
  {
    SurfaceFace face;
    face.tetrahedron = tetra.neighbours[i];
    face.triangle = tetra.triangles[i];
    if (triangles_[face.triangle].is_surface)
      face.curvature = deadFace;
    else
      face.curvature = circumcurvature(tetrahedra_[face.tetrahedron], face.triangle);
    triangles_[face.triangle].surface_face_cached = face;
    surface_.insert(face);
  }
  growSurface(maxCurvature);
}

void ConcaveHull::growOutwards(double maxCurvature)
{
  bool found = false;
  for (auto &tetra : tetrahedra_)
  {
    if (insideTetrahedron(centre_, tetra))
    {
      growOutwards(tetra, maxCurvature);
      found = true;
    }
  }
  if (!found)
    std::cout << "could not find a tetrahedron at the point cloud mean location" << std::endl;
  convertToMesh();
}

// starting with the outer (convex) surface mesh, grow inwards up to the maxCurvature value
void ConcaveHull::growInwards(double maxCurvature)
{
  surface_.clear();
  // find the surface triangles...
  for (int i = 0; i < (int)triangles_.size(); i++)
  {
    Triangle &tri = triangles_[i];
    if (tri.is_surface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra_[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      face.triangle = i;
      face.curvature = circumcurvature(tetrahedra_[face.tetrahedron], face.triangle);
      triangles_[i].surface_face_cached = face;
      surface_.insert(face);
    }
  }
  growSurface(maxCurvature);
  convertToMesh();
}

void ConcaveHull::growInDirection(double maxCurvature, const Eigen::Vector3d &dir)
{
  surface_.clear();
  // find the surface triangles...
  for (int i = 0; i < (int)triangles_.size(); i++)
  {
    Triangle &tri = triangles_[i];
    if (tri.is_surface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra_[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      if (face.tetrahedron < 0)
        std::cout << "bad face tetrahedron" << std::endl;
      Eigen::Vector3d mid(0, 0, 0);
      for (int j = 0; j < 4; j++) mid += vertices_[tetrahedra_[face.tetrahedron].vertices[j]] / 4.0;
      Eigen::Vector3d normal = (vertices_[tri.vertices[2]] - vertices_[tri.vertices[0]])
                                 .cross(vertices_[tri.vertices[1]] - vertices_[tri.vertices[0]]);
      if ((mid - vertices_[tri.vertices[0]]).dot(normal) < 0.0)
        normal = -normal;
      if (normal.dot(dir) < 0.0)  // downwards facing
        continue;
      face.triangle = i;
      for (int j = 0; j < 3; j++)
      {
        vertex_on_surface_[tri.vertices[j]] = true;
        edges_[tri.edges[j]].has_had_face = true;
      }
      face.curvature = circumcurvature(tetrahedra_[face.tetrahedron], face.triangle);
      triangles_[i].surface_face_cached = face;
      surface_.insert(face);
    }
  }
  growSurface(maxCurvature);
  convertToMesh();
}

void ConcaveHull::convertToMesh()
{
  mesh_.vertices() = vertices_;
  int num_bads = 0;
  for (auto &face : surface_)
  {
    Eigen::Vector3d centroid(0, 0, 0);
    ray::ConcaveHull::Tetrahedron &tetra = tetrahedra_[face.tetrahedron];
    Eigen::Vector3i tri_verts = triangles_[face.triangle].vertices;
    if (tri_verts[0] == -1)
      std::cout << "bad vertices in the surface" << std::endl;
    if (tetra.vertices[0] != -1)
    {
      for (int i = 0; i < 4; i++) centroid += vertices_[tetra.vertices[i]] / 4.0;
      Eigen::Vector3d vs[3];
      for (int i = 0; i < 3; i++) vs[i] = vertices_[tri_verts[i]];
      Eigen::Vector3d normal = (vs[2] - vs[0]).cross(vs[1] - vs[0]);
      if ((centroid - vs[0]).dot(normal) < 0.0)
        std::swap(tri_verts[1], tri_verts[2]);
    }
    else
      num_bads++;
    mesh_.indexList().push_back(tri_verts);
  }
  if (num_bads > 0)
    std::cout << "number of surfaces that didn't have enough information to orient: " << num_bads << std::endl;
}
}  // namespace ray
#endif
