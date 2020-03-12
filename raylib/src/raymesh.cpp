// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
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
      if ((contactPoint - corners[i]).dot(side) > 0.0)
        return false;
    }
    return true;
  }
  bool intersectsCube(const Vector3d &cubeMin, double cubeWidth)
  {
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
    for (int j = 0; j<3; j++)
      triangles[i].corners[j] = vertices[indexList[i][j]];
    triangles[i].tested = false;
    triangles[i].normal = (triangles[i].corners[1]-triangles[i].corners[0]).cross(triangles[i].corners[2]-triangles[i].corners[0]).normalized();

    double dir = sgn(offset);
    Vector3d sides[3];
    for (int j = 0; j<3; j++)
      sides[j] = (triangles[i].corners[(j+1)%3] - triangles[i].corners[j]).cross(triangles[i].normal);
    for (int j = 0; j<3; j++)
    {
      Vector3d normal = normals[indexList[i][j]];
      double d1 = normal.dot(sides[j]);
      double d2 = normal.dot(sides[(j+2)%3]);
      
      if (d1*dir < 0.0 && d2*dir < 0.0)
        normal = triangles[i].normal;
      else if (d1*dir < 0.0)
        normal -= sides[j]*d1/sides[j].squaredNorm();
      else if (d2*dir < 0.0)
        normal -= sides[(j+2)%3]*d2/sides[(j+2)%3].squaredNorm();
      normal.normalize();
      triangles[i].corners[j] += normal * offset;  
      boxMin = minVector(boxMin, triangles[i].corners[j]);
      boxMax = maxVector(boxMax, triangles[i].corners[j]);   
    }
  }

#if 0
  Mesh temp;
  for (int i = 0; i<(int)triangles.size(); i++)
  {
    temp.vertices.push_back(triangles[i].corners[0]);
    temp.vertices.push_back(triangles[i].corners[1]);
    temp.vertices.push_back(triangles[i].corners[2]);
    temp.indexList.push_back(Vector3i(3*i,3*i+1,3*i+2));
  }
  writePlyMesh("test.ply", temp);
#endif

  // Thirdly, put the triangles into a grid
  double voxelWidth = 1.0;
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
          if (tri.intersectsCube(boxMin + voxelWidth*Vector3d(x,y,z), voxelWidth))
            grid.cell(x,y,z).data.push_back(&tri);
        }
      }
    }
  }

  // Fourthly, drop each end point downwards to decide whether it is inside or outside..
  for (int r = 0; r<(int)cloud.ends.size(); r++)
  {
    int insides = 0;
    int outsides = 0;
    for (int dir = -1; dir<=1; dir += 2)
    {
      Vector3d start = (cloud.ends[r] - boxMin)/voxelWidth;
      Vector3i index(start[0], start[1], start[2]);
      int endI = dir < 0 ? 1 : grid.dims[2];
      vector<Triangle *> trisTested;
      for (int z = index[2]; (z*dir)<endI; z+=dir)
      {
        auto &tris = grid.cell(index[0], index[1], z).data;
        double minDepth = 10.0;
        Vector3d minNormal;
        for (auto &tri: tris)
        {
          if (tri->tested)
            continue;
          tri->tested = true;
          trisTested.push_back(tri);
          double depth;
          if (tri->intersectsRay(cloud.ends[r], cloud.ends[r] + (double)dir*Vector3d(0.0,0.0,1e4), depth))
          {
            if (depth < minDepth)
            {
              minDepth = depth;
              minNormal = tri->normal;
            }
          }
        }    
        if (minDepth < 10.0)
        {
          if (minNormal[2]*(double)dir < 0.0)
            insides++;
          else
            outsides++;
          break;
        }
      }
      for (auto &tri: trisTested)
        tri->tested = false;
      if (dir==1)
      {
        bool doOutside;
        if (offset >= 0.0)
          doOutside = outsides > 0 || insides == 0;
        else
          doOutside = insides == 0;
        
        if (doOutside)
        {
          outside.starts.push_back(cloud.starts[r]);
          outside.ends.push_back(cloud.ends[r]);
          outside.times.push_back(cloud.times[r]);
        }
        else
        {
          inside.starts.push_back(cloud.starts[r]);
          inside.ends.push_back(cloud.ends[r]);
          inside.times.push_back(cloud.times[r]);
        }     
      }
    }
  }
  cout << "inside: " << inside.starts.size() << ", outside: "  << outside.starts.size() << ". Total " << outside.starts.size() + inside.starts.size()  << " out of " << cloud.ends.size() << " rays" << endl;
#if 0
  // Fifthly, split each ray as it crosses each triangle. We count the number of ins and outs to allow overlapping meshes
  for (int r = 0; r<(int)cloud.ends.size(); r++)
  {
    Vector3d dir = cloud.ends[r] - cloud.starts[r];
    Vector3d dirSign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (cloud.starts[r] - boxMin)/voxelWidth;
    Vector3d end = (cloud.ends[r] - boxMin)/voxelWidth;
    Vector3i index(start[0], start[1], start[2]);
    Vector3i endIndex(end[0], end[1], end[2]);
    for (;;)
    {
      auto &tris = grid.cell(index[0], index[1], index[2]).data;
      for (auto &tri: tris)
        tri.cropRay(cloud.starts[r], cloud.ends[r]);
      if (index == endIndex)
        break;
      Vector3d mid = boxMin + voxelWidth*Vector3d(index[0]+0.5, index[1]+0.5, index[2]+0.5);
      Vector3d nextBoundary = mid + 0.5*voxelWidth*dirSign;
      Vector3d delta = nextBoundary - cloud.starts[r];
      Vector3d d(delta[0]/dir[0], delta[1]/dir[1], delta[2]/dir[2]);
      if (d[0] < d[1] && d[0]<d[2])
        index[0] += dirSign[0];
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += dirSign[1];
      else
        index[2] += dirSign[2];
    }
  }
#endif
}
