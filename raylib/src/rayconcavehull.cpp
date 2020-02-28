#include "rayconcavehull.h"
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

using namespace RAY;

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
  FD_SET (STDIN_FILENO, &rdfs);
 
  select(STDIN_FILENO+1, &rdfs, NULL, NULL, &tv);
  return (bool)FD_ISSET(STDIN_FILENO, &rdfs);
}
#else // windows?
#include <conio.h>
#endif

static const double deadFace = 1e10;
#if defined(DEBUGDRAW)
static DebugDraw *debug = NULL;
#endif

class Hasher
{
public:
  size_t operator() (const Vector2i& key) const
  {
    return key[0] + key[1];
  }
};

ConcaveHull::ConcaveHull(const vector<Vector3d> &points)
{
#if defined(DEBUGDRAW)
  debug = new DebugDraw();
#endif
  centre = mean(points);
  unordered_map<Vector2i, int, Hasher> edgeLookup(points.size()*2);
  
  cout << "number of points: " << points.size() << endl;
  vertices = points;
  vertexOnSurface.resize(vertices.size());
  for (int i = 0; i<(int)vertexOnSurface.size(); i++)
    vertexOnSurface[i] = false;
  
  vector<double> coordinates(points.size() * 3);
  for (int i = 0; i < (int)points.size(); i++) 
  {
    coordinates[3*i + 0] = points[i][0];
    coordinates[3*i + 1] = points[i][1];
    coordinates[3*i + 2] = points[i][2];
  }

  orgQhull::Qhull hull;
  hull.setOutputStream(&cout);
  hull.runQhull("", 3, points.size(), coordinates.data(), "d Qbb Qt");

  orgQhull::QhullFacetList facets = hull.facetList();
  int maxFacets = 0;
  for (const orgQhull::QhullFacet &f : facets) 
    maxFacets = max(maxFacets, f.id()+1);
  cout << "number of total facets: " << facets.size() << endl;
  tetrahedra.resize(maxFacets);
  
  int maxTris = 0; 
  for (const orgQhull::QhullFacet &f : facets) 
  {
    if (f.isUpperDelaunay())
      continue;
    qh_makeridges(hull.qh(), f.getFacetT());
    for (const orgQhull::QhullRidge &r : f.ridges())
      maxTris = max(maxTris, r.id()+1);
  }
  triangles.resize(maxTris);
  cout << "maximum number of triangles: " << maxTris << endl;
  
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
      Vector3d vec(data[0], data[1], data[2]);
      tetra.vertices[i++] = v.point().id();
      if ((vertices[v.point().id()] - vec).squaredNorm() > 1e-8)
        cout << "vertex data doesn't match its id" << endl;
    }

    orgQhull::QhullRidgeSet ridges = f.ridges();
    if (ridges.size() != 4)
      cout << "bad number of ridges: " << ridges.size() << endl;
    i = 0;
    for (const orgQhull::QhullRidge &r : ridges)
    {
      if (r.vertices().size() != 3)
        cout << "bad ridge size: " << r.vertices().size() << endl;
      tetra.neighbours[i] = r.topFacet() == f ? r.bottomFacet().id() : r.topFacet().id();
      tetra.triangles[i++] = r.id();
      if (r.id() > maxTris)
        cout << "bad rid" << endl;
      int j = 0;
      int rid = r.id();
      if (rid >= (int)triangles.size() || r.id() < 0)
        cout << "bag bad" << endl;
      if (triangles[rid].valid())
      {
        if (triangles[rid].tetrahedra[0] != r.topFacet().id() && triangles[rid].tetrahedra[0] != r.bottomFacet().id())
          cout << "replacing with bad data" << endl;
      }
      else
      {
        triangles[rid].tetrahedra[0] = r.topFacet().id();
        triangles[rid].tetrahedra[1] = r.bottomFacet().id();
        triangles[rid].isSurface = r.topFacet().isUpperDelaunay() != r.bottomFacet().isUpperDelaunay();
        if (triangles[rid].tetrahedra[0] != f.id() && triangles[rid].tetrahedra[1] != f.id())
          cout << "bad too" << endl;
        for (const orgQhull::QhullVertex &v : r.vertices()) 
        {
          int vid = v.point().id();
          if (vid != tetra.vertices[0] && vid != tetra.vertices[1] && vid != tetra.vertices[2] && vid != tetra.vertices[3])
            cout << "bad" << endl;
          triangles[rid].vertices[j++] = vid;
        }
        for (int i = 0; i<3; i++)
        {
          int a = triangles[rid].vertices[i];
          int b = triangles[rid].vertices[(i+1)%3];
          Vector2i v(min(a,b), max(a,b));
          const auto &res = edgeLookup.find(v);
          if (res == edgeLookup.end())
          {
            triangles[rid].edges[i] = edges.size();
            edgeLookup.insert({v, edges.size()});
            edges.push_back(Edge(v[0], v[1]));
          }
          else
            triangles[rid].edges[i] = res->second;
        }
      }
    }
    
    if (f.id() >= (int)tetrahedra.size())
      cout << "bad" << endl;
    tetrahedra[f.id()] = tetra;
    
    c++;
  }
  cout << "number of tetrahedrons: " << c << endl;
}

double ConcaveHull::circumcurvature(const ConcaveHull::Tetrahedron &tetra, int triangleID)
{
  Triangle &triangle = triangles[triangleID];
  Vector3d vs[4];
  for (int i = 0; i<4; i++)
    vs[i] = vertices[tetra.vertices[i]];

  Vector3d m1 = (vs[0] + vs[1])*0.5;
  Vector3d m2 = (vs[0] + vs[2])*0.5;
  Vector3d normal = (vs[2]-vs[0]).cross(vs[1] - vs[0]);
  Vector3d up = (vs[2] - vs[0]).cross(normal);
  
  // defines the line midway between vs[0], vs[1] and vs[2]
  double t1 = (m1 - m2).dot(vs[1]-vs[0])/up.dot(vs[1]-vs[0]);
  Vector3d midbase = m2 + up * t1;
  
  // defines the plane midway between vs[0] and vs[3]
  Vector3d m = (vs[3] + vs[0])*0.5;
  Vector3d dir = vs[3] - vs[0];
  
  // intersection of line and plane is midway between all four points
  double t = (m-midbase).dot(dir) / normal.dot(dir);
  Vector3d circumcentre = midbase + normal * t;
  Vector3d centre = (vs[0] + vs[1] + vs[2] + vs[3])*0.25;
  
  // return the distance from a vertex to this circumcentre
  double circumradius = (circumcentre - vs[0]).norm();
  
  Vector3d cs[3];
  for (int i = 0; i<3; i++)
    cs[i] = vertices[triangle.vertices[i]];

  Vector3d triNormal = (cs[2]-cs[0]).cross(cs[1] - cs[0]);
  double circumcentreSide = (circumcentre - cs[0]).dot(triNormal);
  if (circumcentreSide * (centre - cs[0]).dot(triNormal) > 0.0) 
  {
    // we are ready to make this a dead face, but before we do, it is possible that the tetrahedron in question includes an existing surface face
    // in which case, the circumcentre might not in fact be above the current surface.
    
    int faceIntersects = -1;
    int numFaceIntersects = 0;
    for (int j = 0; j<4; j++)
    {
      if (tetra.triangles[j] == triangleID)
        continue;
      if (triangles[tetra.triangles[j]].surfaceFaceCached.triangle != -1)
      {
        numFaceIntersects++;
        faceIntersects = j;      
      }
    }
    if (numFaceIntersects == 1)
    {
      int otherFace = tetra.triangles[faceIntersects];
      Triangle tri = triangles[otherFace];
      Vector3d ds[3];
      for (int i = 0; i<3; i++)
        ds[i] = vertices[tri.vertices[i]];

      Vector3d triNormal2 = (ds[2]-ds[0]).cross(ds[1] - ds[0]);
      double circumcentreSide2 = (circumcentre - ds[0]).dot(triNormal2);
      bool aboveOtherFace = circumcentreSide2 * (centre - ds[0]).dot(triNormal2) > 0.0;
      if (aboveOtherFace)
        return deadFace;
    }
    else 
      return deadFace;
  }

  return 1.0/circumradius;
}

static int newTriCount = 0;

bool ConcaveHull::growFront(double maxCurvature)
{
  SurfaceFace face = *surface.begin();
  if (face.curvature == deadFace || face.curvature > maxCurvature)
    return false;
  int vertexI = 0;
  for (int i = 0; i<4; i++)
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
  bool intersects = vertexOnSurface[newVertex];
  int faceIntersects = -1;
  int numFaceIntersects = 0;
  const SurfaceFace *faceIntersectTri = NULL;
  for (int j = 0; j<4; j++)
  {
    if (tetra.triangles[j] == face.triangle)
      continue;
    if (triangles[tetra.triangles[j]].surfaceFaceCached.triangle != -1)
    {
      numFaceIntersects++;
      faceIntersects = j;      
      faceIntersectTri = &triangles[tetra.triangles[j]].surfaceFaceCached;
    }
  }
 
  surface.erase(surface.begin());
  if (numFaceIntersects == 1)
  {
    intersects = false;
    int otherVertex = 0;
    for (int i = 0; i<3; i++)
    {
      int v = triangles[face.triangle].vertices[i];
      Triangle &tri = triangles[tetra.triangles[faceIntersects]];
      if (v != tri.vertices[0] && v != tri.vertices[1] && v != tri.vertices[2])
        otherVertex = v;
    }
    int tri2 = tetra.triangles[(faceIntersects+1)%3];
    if (tri2 == face.triangle)
      tri2 = tetra.triangles[(faceIntersects+2)%3];
    int v0 = min(otherVertex, newVertex);
    int v1 = max(otherVertex, newVertex);
    for (int i = 0; i<3; i++)
    {
      Edge &edge = edges[triangles[tri2].edges[i]];      
      if (edge.vertices[0] == v0 && edge.vertices[1] == v1 && edge.hasHadFace)
        intersects = true;
    } 
    if (!intersects)
      surface.erase(*faceIntersectTri);
  }
  if (intersects)
  {
    face.curvature = deadFace;
    triangles[face.triangle].surfaceFaceCached = face;
    surface.insert(face); // put it at the back of the queue
    return true;
  }

  newTriCount++;
  for (int i = 0; i<4; i++)
  {
    if (tetra.triangles[i] == face.triangle || (numFaceIntersects == 1 && i==faceIntersects))
      continue;
    SurfaceFace newFace;
    
    newFace.tetrahedron = tetra.neighbours[i];
    newFace.triangle = tetra.triangles[i];
    for (int j = 0; j<3; j++)
      vertexOnSurface[triangles[newFace.triangle].vertices[j]] = true;
    double grad;
    if (triangles[newFace.triangle].isSurface)
      newFace.curvature = grad = deadFace;
    else
      newFace.curvature = circumcurvature(tetrahedra[newFace.tetrahedron], newFace.triangle);
    triangles[newFace.triangle].surfaceFaceCached = newFace;
    for (int j = 0; j<3; j++)
      edges[triangles[newFace.triangle].edges[j]].hasHadFace = true;
    surface.insert(newFace);
  }
  return true;
}

void ConcaveHull::growSurface(double maxCurvature)
{
  do
  {
    if (!(newTriCount%1600))
    {
      cout << "max curvature of structure: " << surface.begin()->curvature << endl;
      vector< vector<Vector3d> > tris;
      for (auto &tri: surface)
      {
        vector<Vector3d> corners(3);
        for (int i = 0; i<3; i++)
          corners[i] = vertices[triangles[tri.triangle].vertices[i]];
        tris.push_back(corners);
      }
#if defined(DEBUGDRAW)
      debug->drawTriangles(tris, 0.5);
#endif      
      if (kbhit())
        return;
    }
  } while (growFront(maxCurvature));
}

// starting with given tetrahedron, grow it outwards to achieve a maximum curvature
void ConcaveHull::growOutwards(const ConcaveHull::Tetrahedron &tetra, double maxCurvature)
{
  surface.clear();
  for (int i = 0; i<4; i++)
  {
    SurfaceFace face;
    face.tetrahedron = tetra.neighbours[i];
    face.triangle = tetra.triangles[i];
    if (triangles[face.triangle].isSurface)
      face.curvature = deadFace;
    else
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
    triangles[face.triangle].surfaceFaceCached = face;
    surface.insert(face);
  }
  growSurface(maxCurvature);
}

void ConcaveHull::growOutwards(double maxCurvature)
{
  bool found = false;
  for (auto &tetra: tetrahedra)
  {
    if (insideTetrahedron(centre, tetra))
    {
      growOutwards(tetra, maxCurvature);
      found = true;
    }
  }
  if (!found)
    cout << "could not find a tetrahedron at the point cloud mean location" << endl;
}

// starting with the outer (convex) surface mesh, grow inwards up to the maxCurvature value
void ConcaveHull::growInwards(double maxCurvature)
{
  surface.clear();
  // find the surface triangles...
  for (int i = 0; i<(int)triangles.size(); i++)
  {
    Triangle &tri = triangles[i];
    if (tri.isSurface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      face.triangle = i;
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
      triangles[i].surfaceFaceCached = face;
      surface.insert(face);
    }
  }
  growSurface(maxCurvature);
}

void ConcaveHull::growInDirection(double maxCurvature, const Vector3d &dir)
{
  surface.clear();
  // find the surface triangles...
  for (int i = 0; i<(int)triangles.size(); i++)
  {
    Triangle &tri = triangles[i];
    if (tri.isSurface)
    {
      SurfaceFace face;
      face.tetrahedron = tetrahedra[tri.tetrahedra[0]].valid() ? tri.tetrahedra[0] : tri.tetrahedra[1];
      if (face.tetrahedron < 0)
        cout << "bad face tetrahedron" << endl;
      Vector3d mid(0,0,0);
      for (int j = 0; j<4; j++)
        mid += vertices[tetrahedra[face.tetrahedron].vertices[j]]/4.0;
      Vector3d normal = (vertices[tri.vertices[2]] - vertices[tri.vertices[0]]).cross(vertices[tri.vertices[1]] - vertices[tri.vertices[0]]);
      if ((mid - vertices[tri.vertices[0]]).dot(normal) < 0.0)
        normal = -normal;
      if (normal.dot(dir) < 0.0) // downwards facing
        continue;
      face.triangle = i;
      for (int j = 0; j<3; j++)
      {
        vertexOnSurface[tri.vertices[j]] = true;
        edges[tri.edges[j]].hasHadFace = true;
      }
      face.curvature = circumcurvature(tetrahedra[face.tetrahedron], face.triangle);
      triangles[i].surfaceFaceCached = face;
      surface.insert(face);
    }
  }
  growSurface(maxCurvature);
}

