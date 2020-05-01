// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"

#include "raylaz.h"
#include "rayply.h"
#include "raytrajectory.h"
#include "rayunused.h"

#include <set>

using namespace std;
using namespace Eigen;
using namespace RAY;

struct Triangle
{
  Vector3d corners[3];
  Vector3d normal;
  bool tested;
  bool intersectsRay(const Vector3d &rayStart, const Vector3d &rayEnd, double &depth)
  {
    // 1. plane test:
    double d1 = (rayStart - corners[0]).dot(normal);
    double d2 = (rayEnd - corners[0]).dot(normal);
    if (d1*d2 > 0.0)
      return false;

    depth = d1/(d1-d2);
    Vector3d contactPoint = rayStart + (rayEnd - rayStart) * depth;

    // next we have to test every sideways direction
    for (int i = 0; i<3; i++)
    {
      Vector3d side = (corners[(i+1)%3] - corners[i]).cross(normal);
      if ((contactPoint - corners[i]).dot(side) >= 0.0)
        return false;
    }
    return true;
  }
  double distSqrToPoint(const Vector3d &point)
  {
    Vector3d pos = point - normal*(point - corners[0]).dot(normal);
    bool outs[3];
    double ds[3];
    Vector3d sides[3];
    for (int i = 0; i<3; i++)
    {
      sides[i] = (corners[(i+1)%3] - corners[i]).cross(normal);
      ds[i] = (pos - corners[i]).dot(sides[i]);
      outs[i] = ds[i]>0.0;
    }
    if (outs[0] && outs[1])
      pos = corners[1];
    else if (outs[1] && outs[2])
      pos = corners[2];
    else if (outs[2] && outs[0])
      pos = corners[0];
    else if (outs[0])
      pos -= sides[0]*ds[0]/sides[0].squaredNorm();
    else if (outs[1])
      pos -= sides[1]*ds[1]/sides[1].squaredNorm();
    else if (outs[2])
      pos -= sides[2]*ds[2]/sides[2].squaredNorm();
    return (point - pos).squaredNorm();
  }
  bool intersectsCube(const Vector3d &cubeMin, double cubeWidth)
  {
    RAYLIB_UNUSED(cubeMin);
    RAYLIB_UNUSED(cubeWidth);
    // TODO: fill in
    return true;
  }
};

void Mesh::splitCloud(const Cloud &cloud, double offset, Cloud &inside, Cloud &outside)
{
  // Firstly, find the average vertex normals
  vector<Vector3d> normals(vertices.size());
  for (auto &normal: normals)
    normal.setZero();
  for (auto &index: indexList)
  {
    Vector3d normal = (vertices[index[1]]-vertices[index[0]]).cross(vertices[index[2]]-vertices[index[0]]);
    for (int i = 0; i<3; i++)
      normals[index[i]] += normal;
  }
  for (auto &normal: normals)
    normal.normalize();

  // convert to separate triangles for convenience
  vector<Triangle> triangles(indexList.size());
  Vector3d boxMin(1e10,1e10,1e10), boxMax(-1e10,-1e10,-1e10);
  for (int i = 0; i<(int)indexList.size(); i++)
  {
    Triangle &tri = triangles[i];
    for (int j = 0; j<3; j++)
      tri.corners[j] = vertices[indexList[i][j]];
    tri.tested = false;
    tri.normal = (tri.corners[1]-tri.corners[0]).cross(tri.corners[2]-tri.corners[0]).normalized();
    for (int j = 0; j<3; j++)
    {
      boxMin = minVector(boxMin, tri.corners[j]);
      boxMax = maxVector(boxMax, tri.corners[j]);   
    }
  }

  // Thirdly, put the triangles into a grid
  double voxelWidth = 1.0;
  vector<int> insideIndices;
  {
    int insideVal = offset >= 0.0 ? 1 : 0;
    Grid<Triangle *> grid(boxMin, boxMax, voxelWidth);
    for (auto &tri: triangles)
    {
      Vector3d triMin = (minVector(tri.corners[0], minVector(tri.corners[1], tri.corners[2])) - boxMin)/voxelWidth;
      Vector3d triMax = (maxVector(tri.corners[0], maxVector(tri.corners[1], tri.corners[2])) - boxMin)/voxelWidth;
      for (int x = (int)triMin[0]; x<=(int)triMax[0]; x++)
      {
        for (int y = (int)triMin[1]; y<=(int)triMax[1]; y++)
        {
          for (int z = (int)triMin[2]; z<=(int)triMax[2]; z++)
          {
 //           if (tri.intersectsCube(boxMin + voxelWidth*Vector3d(x,y,z), voxelWidth))
              grid.insert(x,y,z, &tri);
          }
        }
      }
    }

    // Fourthly, drop each end point downwards to decide whether it is inside or outside..
    vector<Triangle *> trisTested;
    trisTested.reserve(100.0);
    for (int r = 0; r<(int)cloud.ends.size(); r++)
    {
      int intersections = 0;
     // for (int dir = -1; dir<=1; dir += 2)
      int dir = -1;
      {
        Vector3d start = (cloud.ends[r] - boxMin)/voxelWidth;
        Vector3i index(start.cast<int>());
        int endI = dir < 0 ? 0 : grid.dims[2]-1;
        trisTested.clear();
        for (int z = clamped(index[2], 0, grid.dims[2]-1); (z*dir)<=endI; z+=dir)
        {
          auto &tris = grid.cell(index[0], index[1], z).data;
          for (auto &tri: tris)
          {
            if (tri->tested)
              continue;
            tri->tested = true;
            trisTested.push_back(tri);
            double depth;
            if (tri->intersectsRay(cloud.ends[r], cloud.ends[r] + (double)dir*Vector3d(0.0,0.0,1e3), depth))
              intersections++;
          }    
        }
        for (auto &tri: trisTested)
          tri->tested = false;    
      }
      if ((intersections%2) == insideVal) // inside
        insideIndices.push_back(r);
    }
  }
  cout << insideIndices.size()<<"/"<<cloud.ends.size() << " inside mesh" << endl;

  // what if offset is negative?
  // we have to ...
  if (offset != 0.0)
  {
    // Thirdly, put the triangles into a grid
    double voxelWidth = 1.0;
    Grid<Triangle *> grid2(boxMin, boxMax, voxelWidth);
    for (int i = 0; i<(int)indexList.size(); i++)
    {
      if (!(i%100000))
        cout << "filling volumes " << i << "/" << indexList.size() << endl;
      Triangle &tri = triangles[i];
      Vector3d extrudedCorners[3];
      for (int j = 0; j<3; j++)
        extrudedCorners[j] = tri.corners[j] + normals[indexList[i][j]]*offset;
      
      Vector3d triMin = minVector(tri.corners[0], minVector(tri.corners[1], tri.corners[2]));
      Vector3d triMax = maxVector(tri.corners[0], maxVector(tri.corners[1], tri.corners[2]));
      Vector3d triMin2 = minVector(extrudedCorners[0], minVector(extrudedCorners[1], extrudedCorners[2]));
      Vector3d triMax2 = maxVector(extrudedCorners[0], maxVector(extrudedCorners[1], extrudedCorners[2]));

      triMin = (minVector(triMin, triMin2) - boxMin)/voxelWidth;
      triMax = (maxVector(triMax, triMax2) - boxMin)/voxelWidth;
      for (int x = (int)triMin[0]; x<=(int)triMax[0]; x++)
        for (int y = (int)triMin[1]; y<=(int)triMax[1]; y++)
          for (int z = (int)triMin[2]; z<=(int)triMax[2]; z++)
            grid2.insert(x,y,z, &tri);
    }  
    // now go through the remaining inside points
    vector<int> newInsides;
    double offsetSqr = sqr(offset);
    int p = 0;
    for (auto &r: insideIndices)
    {
      if (!(p++%1000000))
        cout << "checking points " << p-1 << "/" << insideIndices.size() << endl;
      Vector3d pos = (cloud.ends[r] - boxMin)/voxelWidth;
      Vector3i index(pos.cast<int>());
      auto &tris = grid2.cell(index[0], index[1], index[2]).data;
      bool inTri = false;
      for (auto &tri: tris)
      {
        if (tri->distSqrToPoint(cloud.ends[r]) < offsetSqr)
        {
          inTri = true;
          break;
        };
      }
      if (!inTri)
        newInsides.push_back(r);
    }
    cout << "new inside count: " << newInsides.size() << "/" << cloud.ends.size() << endl;
    insideIndices = newInsides;
  }

  vector<bool> insideI(cloud.ends.size());
  int ins = offset >= 0.0;
  for (int i = 0; i<(int)cloud.ends.size(); i++)
    insideI[i] = !ins;
  for (auto &ind: insideIndices)
    insideI[ind] = ins;
  for (int i = 0; i<(int)cloud.ends.size(); i++)
  {
    Cloud &out = insideI[i] ? inside : outside;
    out.starts.push_back(cloud.starts[i]);
    out.ends.push_back(cloud.ends[i]);
    out.times.push_back(cloud.times[i]);
    out.colours.push_back(cloud.colours[i]);
  }
}
