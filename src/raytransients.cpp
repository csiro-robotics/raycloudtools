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
  cout << "raytransients raycloud 0.1 m 3 s - splits out points intersected within 0.1m more than 3 seconds before or after" << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 6)
    usage();

  if (string(argv[3]) != "m" || string(argv[5]) != "s")
    usage();

  string file = argv[1];
  Cloud cloud;
  cloud.load(file);

  double radius = stod(argv[2]);
  double timeDelta = stod(argv[4]);

  Cloud transient;
  Cloud fixed;
  cloud.findTransients(transient, fixed, radius, timeDelta);

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);

  transient.save(fileStub + "_transient.ply");
  fixed.save(fileStub + "_fixed.ply");
  return true;
}
