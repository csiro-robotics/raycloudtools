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
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 5)
    usage();
  
  string file = argv[1];
  Pose pose;
  pose.position = Vector3d(stod(argv[2]), stod(argv[3]), stod(argv[4]));
  pose.rotation = Quaterniond(1.0,0.0,0.0,0.0);

  Cloud cloud;
  cloud.load(file);
  cloud.transform(pose, 0.0);
  cloud.save(file);
  return true;
}
