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
  cout << "Splits a raycloud into the transient rays and the fixed part" << endl;
  cout << "usage:" << endl;
  cout << "raytransients min raycloud 3 s - splits out transient points more than 3 seconds apart from the crossing rays" << endl;
  cout << "              max    - finds negative transients, such as a hallway exposed when a door opens." << endl;
  cout << "              oldest - keeps the oldest geometry when there is a difference over time." << endl;
  cout << "              newest - uses the newest geometry when there is a difference over time." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 5)
    usage();

  if (string(argv[4]) != "s")
    usage();
  string mergeType = argv[1];
  if (mergeType != "min" && mergeType != "max" && mergeType != "oldest" && mergeType != "newest")
    usage();
  string file = argv[2];
  Cloud cloud;
  cloud.load(file);

  double timeDelta = stod(argv[3]);

  Cloud transient;
  Cloud fixed;
  cloud.findTransients(transient, fixed, timeDelta, mergeType);

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);

  transient.save(fileStub + "_transient.ply");
  fixed.save(fileStub + "_fixed.ply");
  return true;
}
