// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raymesh.h"
#include "rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Split a ray cloud relative to the supplied triangle mesh, generating two cropped ray clouds" << endl;
  cout << "usage:" << endl;
  cout << "raysplit raycloud meshfile 0.2  - splits raycloud along mesh offset by 0.2 m along the vertex normals" << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4)
    usage();
  
  string file = argv[1];
  Cloud cloud;
  cloud.load(file);

  string meshFile = argv[2];
  Mesh mesh;
  readPlyMesh(meshFile, mesh);

  double offset = stod(argv[3]);
  Cloud inside, outside;
  mesh.splitCloud(cloud, offset, inside, outside);

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);

  inside.save(fileStub + "_inside.ply");
  outside.save(fileStub + "_outside.ply");
  return true;
}
