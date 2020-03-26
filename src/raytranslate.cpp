// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;


void usage(bool error=false)
{
  cout << "Translate a raycloud" << endl;
  cout << "usage:" << endl;
  cout << "raytranslate raycloud 0 0 1 - translation (x,y,z) in metres" << endl;
  cout << "                      0 0 1 24.3 - optional 4th component translates time" << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 5 && argc != 6)
    usage();
  
  string file = argv[1];
  Pose pose;
  pose.position = Vector3d(stod(argv[2]), stod(argv[3]), stod(argv[4]));
  double timeDelta = 0.0;
  if (argc == 6)
    timeDelta = stod(argv[5]);
  pose.rotation = Quaterniond(1.0,0.0,0.0,0.0);

  Cloud cloud;
  cloud.load(file);
  #if 0 // test bending
  Vector3d bend(pose.position[0], pose.position[1], 0.0);
  Vector3d side = Vector3d(bend[1], -bend[0], 0).normalized();
  for (auto &end: cloud.ends)
    end += bend*sqr(end.dot(side));
  #else
  cloud.transform(pose, timeDelta);
  #endif
  cloud.save(file);
  return true;
}
