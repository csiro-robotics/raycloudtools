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
  cout << "Decimate a ray cloud spatially or temporally" << endl;
  cout << "usage:" << endl;
  cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << endl;
  cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4)
    usage();

  string file = argv[1];
  Cloud cloud;
  cloud.load(file);

  Cloud new_cloud;
  string type = argv[3];
  if (type == "cm")
  {
    cloud.decimate(0.01 * stod(argv[2]));
  }
  else if (type == "rays")
  {
    int decimation = stoi(argv[2]);
    for (int i = 0; i < (int)cloud.ends.size(); i += decimation)
      new_cloud.addRay(cloud, i);
    cloud = new_cloud;
  }
  else
    usage(false);
  string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);
  cloud.save(file_stub + "_decimated.ply");
  return true;
}
