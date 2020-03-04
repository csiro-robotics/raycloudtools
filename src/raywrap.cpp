/* (c) Copyright CSIRO 2013. Author: Thomas Lowe
   This software is provided under the terms of Schedule 1 of the license agreement between CSIRO, 3DLM and GeoSLAM.
*/
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <getopt.h>

#include "rayutils.h"
#include "raycloud.h"
#include "rayconcavehull.h"
#include "rayconvexhull.h"
#include "rayply.h"
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Extracts the ground surface as a mesh." << endl;
  cout << "usage:" << endl;
  cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards" << endl;
  cout << "                               the 1.0 is the maximum curvature to bend to" << endl;
  cout << "--full                       - the full (slower) method accounts for overhangs." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc < 4 || argc>5)
    usage();

#if defined(DEBUGDRAW)
  ros::init(argc, argv, "ConcaveHull");
#endif
  string file = argv[1];
  string type = argv[2];
  double maximumCurvature = stod(argv[3]);
  bool overhangs = false;
  if (argc == 5)
  {
    if (string(argv[4])=="--full" or string(argv[4])=="-f")
      overhangs = true;
    else
      usage();    
  }
    
  Cloud cloud;
  cloud.load(file);
  if (file.substr(file.length()-4)==".ply")
    file = file.substr(0,file.length()-4);

  if (overhangs)
  {
    ConcaveHull concaveHull(cloud.ends);
    if (type == "inwards")
      concaveHull.growInwards(maximumCurvature);
    else if (type == "outwards")
      concaveHull.growOutwards(maximumCurvature);
    else if (type == "upwards")
      concaveHull.growUpwards(maximumCurvature);
    else if (type == "downwards")
      concaveHull.growTopDown(maximumCurvature);
    else
      usage();
    
    
    vector<Vector3i> tris;
    int numBads = 0;
    for (auto &face: concaveHull.surface)
    {
      Vector3d centroid(0,0,0);
      ConcaveHull::Tetrahedron &tetra = concaveHull.tetrahedra[face.tetrahedron];
      Vector3i triVerts = concaveHull.triangles[face.triangle].vertices;
      if (triVerts[0] == -1)
        cout << "bad vertices in the surface" << endl;
      if (tetra.vertices[0] != -1)
      {
        for (int i = 0; i<4; i++)
          centroid += concaveHull.vertices[tetra.vertices[i]] / 4.0;
        Vector3d vs[3];
        for (int i = 0; i<3; i++)
          vs[i] = concaveHull.vertices[triVerts[i]];
        Vector3d normal = (vs[2]-vs[0]).cross(vs[1]-vs[0]);
        if ((centroid - vs[0]).dot(normal) < 0.0)
          swap(triVerts[1], triVerts[2]);
      }
      else
        numBads++;
      tris.push_back(triVerts);
    }
    if (numBads > 0)
      cout << "number of surfaces that didn't have enough information to orient: " << numBads << endl;
    writePlyMesh(file + "_mesh.ply", concaveHull.vertices, tris, true);
  }
  else
  {
    ConvexHull convexHull(cloud.ends);
    if (type == "inwards")
      convexHull.growInwards(maximumCurvature);
    else if (type == "outwards")
      convexHull.growOutwards(maximumCurvature);
    else if (type == "upwards")
      convexHull.growUpwards(maximumCurvature);
    else if (type == "downwards")
      convexHull.growTopDown(maximumCurvature);
    else 
      usage(); 

    vector<Vector3i> tris(convexHull.triangles.size());
    vector<Vector3d> vertices;
    for (int i = 0; i<(int)convexHull.triangles.size(); i++)
    {
      tris[i] = Vector3i(3*i, 3*i + 1, 3*i + 2); // TODO: currently the triangles are indexed independently, I need to improve this
      vertices.push_back(convexHull.triangles[i].vertices[0]);
      vertices.push_back(convexHull.triangles[i].vertices[1]);
      vertices.push_back(convexHull.triangles[i].vertices[2]);
    }
    writePlyMesh(file + "_mesh.ply", vertices, tris, true);     
  }
  

  cout << "Completed, output: " << file << "_mesh.ply" << endl;
}


