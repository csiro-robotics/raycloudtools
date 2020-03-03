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
  cout << "Transform a raycloud as a translation, rotation and time shift" << endl;
  cout << "usage:" << endl;
  cout << "raytransform raycloud 0 0 1  30 0 0  13.5 - translation (x,y,z) rotation (rx,ry,rz) time shift (t):" << endl;
  cout << "                                            so this example raises the cloud by 1 metre," << endl;
  cout << "                                            rotates 30 degrees around the (1,0,0) x-axis," << endl;
  cout << "                                            and shifts the time forwards 13.5 seconds." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 9)
    usage();
  
  string file = argv[1];
  Pose pose;
  pose.position = Vector3d(stod(argv[2]), stod(argv[3]), stod(argv[4]));
  Vector3d axis(stod(argv[5]), stod(argv[6]), stod(argv[7]));
  double angle = axis.norm();
  axis /= angle;
  pose.rotation = Quaterniond(AngleAxisd(angle, axis));
  double timeDelta = stod(argv[8]);

  Cloud cloud;
  cloud.load(file);
  cloud.transform(pose, timeDelta);
  cloud.save(file);
  return true;
}
