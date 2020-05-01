// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(int exitCode = 0)
{
  cout << "Split a ray cloud relative to the supplied triangle mesh, generating two cropped ray clouds" << endl;
  cout << "usage:" << endl;
  cout << "raysplit raycloud pos 10,0,0             - splits along x axis" << endl;
  cout << "                  time 10000             - splits at given acquisition time" << endl;
  cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << endl;
  cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << endl;
  cout << "                  meshfile distance 0.2  - splits raycloud at 0.2m from the meshfile surface" << endl;
  cout << "                  startpos 1,2,3         - splits based on start position, around plane 1,2,3" << endl;
  cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays" << endl;
  cout << "                  range 10               - splits out rays more than 10 m long" << endl;
  cout << "                  speed 1.0              - splits out rays when sensor moving above the given speed" << endl;
  exit(exitCode);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 5)
    usage();
  
  string file = argv[1];
  Cloud cloud;
  cloud.load(file);

  Cloud inside, outside;
  if (argc == 5)
  {
    string meshFile = argv[2];
    Mesh mesh;
    readPlyMesh(meshFile, mesh);

    double offset = stod(argv[4]);
    mesh.splitCloud(cloud, offset, inside, outside);
  }
  else
  {
    string parameter = string(argv[2]);
    if (parameter == "time")
    {
      double val = stod(argv[3]);
      cloud.split(inside, outside, 
        [&](int i) -> bool 
        {
          return cloud.times[i] > val; 
        });
    }
    else if (parameter == "alpha")
    {
      double val = stod(argv[3]);
      if (!(val>=0.0 && val <= 1.0))
        usage();
      uint8_t c = uint8_t(255.0*val);
      cloud.split(inside, outside, 
        [&](int i)
        {
          return cloud.colours[i].alpha > c; 
        });
    }
    else if (parameter == "pos")
    {
      stringstream ss(argv[3]);
      Vector3d vec;
      ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, 
        [&](int i)
        {
          return cloud.ends[i].dot(vec) > 1.0; 
        });
    }
    else if (parameter == "startpos")
    {
      stringstream ss(argv[3]);
      Vector3d vec;
      ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, 
        [&](int i)
        {
          return cloud.starts[i].dot(vec) > 0.0; 
        });
    }
    else if (parameter == "raydir")
    {
      stringstream ss(argv[3]);
      Vector3d vec;
      ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, 
        [&](int i)
        {
          Vector3d rayDir = (cloud.ends[i] - cloud.starts[i]).normalized();
          return rayDir.dot(vec) > 0.0; 
        });
    }
    else if (parameter == "colour")
    {
      stringstream ss(argv[3]);
      Vector3d vec;
      ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, 
        [&](int i)
        {
          Vector3d col((double)cloud.colours[i].red/255.0, (double)cloud.colours[i].green/255.0, (double)cloud.colours[i].blue/255.0);
          return col.dot(vec) > 0.0; 
        });
    }
    else if (parameter == "range")
    {
      double val = stod(argv[3]);
      cloud.split(inside, outside, 
        [&](int i)
        {
          return (cloud.starts[i]-cloud.ends[i]).norm() > val; 
        });
    }
    else if (parameter == "speed")
    {
      double val = stod(argv[3]);
      cloud.split(inside, outside, 
        [&](int i)
        {
          if (i==0)
            return false;
          return (cloud.starts[i]-cloud.starts[i-1]).norm()/(cloud.times[i]-cloud.times[i-1]) > val; 
        });
    }
  }

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);

  inside.save(fileStub + "_inside.ply");
  outside.save(fileStub + "_outside.ply");
  return true;
}
