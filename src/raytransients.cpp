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
  cout << "raytransients raycloud 3 s - splits out transient points more than 3 seconds apart from the crossing rays" << endl;
  cout << "              --intersect  - intersection rather than union, this removes negative transients, like " << endl;
  cout << "                             hallways that are temporarily exposed when a door opens. " << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 5)
    usage();

  if (string(argv[3]) != "s")
    usage();
  bool maximal = false;
  if (argc == 5)
  {
    if (string(argv[4]) != "--intersect" && string(argv[4]) != "-i")
      usage();
    maximal = true;
  }

  string file = argv[1];
  Cloud cloud;
  cloud.load(file);

  double timeDelta = stod(argv[2]);

  Cloud transient;
  Cloud fixed;
  cloud.findTransients(transient, fixed, timeDelta, maximal);

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);

  transient.save(fileStub + "_transient.ply");
  fixed.save(fileStub + "_fixed.ply");
  return true;
}
