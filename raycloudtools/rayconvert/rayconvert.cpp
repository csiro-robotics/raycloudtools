// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace ray;

void usage(int exit_code = 0)
{
  cout << "Convert a point cloud and trajectory file into a ray cloud" << endl;
  cout << "usage:" << endl;
  cout << "rayconvert pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << endl;
  cout << "                                            trajectoryfile is a text file in time,x,y,z,qw,qx,qy,qz format"
       << endl;
  cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();

  string point_cloud_file = argv[1];
  string trajectory_file = argv[2];
  Cloud cloud;
  cloud.load(point_cloud_file, trajectory_file);

  string save_file = point_cloud_file.substr(0, point_cloud_file.size() - 4);
  if (point_cloud_file.substr(point_cloud_file.size() - 4) == ".ply")
    save_file += "_raycloud";
  cloud.save(save_file + ".ply");

  return true;
}
