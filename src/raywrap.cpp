/* (c) Copyright CSIRO 2013. Author: Thomas Lowe
   This software is provided under the terms of Schedule 1 of the license agreement between CSIRO, 3DLM and GeoSLAM.
*/
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <getopt.h>

#include "rayutils.h"
#include "rayconcavehull.h"
#include "rayply.h"

void usage(bool error=false)
{
  cout << "Extracts the ground surface as a mesh." << endl;
  cout << "usage:" << endl;
  cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards" << endl;
  cout << "                               the 1.0 is the maximum curvature to bend to" << endl;
  cout << "--fast                       - a fast method that doesn't account for overhangs." << endl;
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
  bool overhangs = true;
  if (argc == 5 && string(argv[4])=="--fast" or string(argv[4])=="-f")
    overhangs = false;
    
  Cloud cloud;
  cloud.load(file);

  if (overhangs)
  {
    ConcaveHull concaveHull(cloud.ends);
    if (type == "inwards")
      concaveHull.growInwards(maximumCurvature);
    else if (type == "outwards")
      concaveHull.growOutwards(maximumCurvature);
    else if (type == "upwards")
      concaveHull.growBottomUp(maximumCurvature);
    else if (type == "downwards")
      concaveHull.growTopDown(maximumCurvature);
    
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
      convexHull.growBottomUp(maximumCurvature);
    else if (type == "downwards")
      convexHull.growTopDown(maximumCurvature); 
    vector<Vector3i> tris(convexHull.vertices.size());
    int c = 0;
    for (int i = 0; i<(int)convexHull.vertices.size(); i++)
      tris[i] = Vector3i(c++,c++,c++); // TODO: currently the triangles are indexed independently, I need to improve this
    writePlyMesh(file + "_mesh.ply", convexHull.vertices, tris, true);     
  }
  

  cout << "Completed, output: " << file << "_mesh.ply" << endl;
}


