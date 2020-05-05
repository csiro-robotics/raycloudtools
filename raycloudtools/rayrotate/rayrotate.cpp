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
using namespace RAY;


void usage(int exitCode = 0)
{
  cout << "Rotate a raycloud about the origin" << endl;
  cout << "usage:" << endl;
  cout << "rayrotate raycloud 30,0,0  - rotation (rx,ry,rz) is a rotation vector in degrees:" << endl;
  cout << "                             so this example rotates the cloud by 30 degrees in the x axis." << endl;
  exit(exitCode);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();
  
  string file = argv[1];
  Pose pose;
  pose.position = Vector3d(0,0,0);

  stringstream ss(argv[2]);
  Vector3d axis;
  ss >> axis[0]; ss.ignore(1); ss>>axis[1]; ss.ignore(1); ss>>axis[2];

  double angle = axis.norm();
  axis /= angle;
  pose.rotation = Quaterniond(AngleAxisd(angle * pi/180.0, axis));

  Cloud cloud;
  cloud.load(file);
  cloud.transform(pose, 0.0);
  cloud.save(file);
  return true;
}
