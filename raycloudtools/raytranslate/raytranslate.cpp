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
  cout << "Translate a raycloud" << endl;
  cout << "usage:" << endl;
  cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << endl;
  cout << "                      0,0,1,24.3 - optional 4th component translates time" << endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();
  
  string file = argv[1];
  Pose pose;
  stringstream ss(argv[2]);
  Vector3d vec;
  ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
  pose.position = vec;

  double time_delta = 0.0;
  if (ss.peek() == ',')
  {
    ss.ignore(1);
    ss >> time_delta;
  }
  pose.rotation = Quaterniond(1.0,0.0,0.0,0.0);

  Cloud cloud;
  cloud.load(file);
  #if 0 // test bending
  Vector3d bend(pose.position[0], pose.position[1], 0.0);
  Vector3d side = Vector3d(bend[1], -bend[0], 0).normalized();
  for (auto &end: cloud.ends)
    end += bend*sqr(end.dot(side));
  #else
  cloud.transform(pose, time_delta);
  #endif
  cloud.save(file);
  return true;
}
