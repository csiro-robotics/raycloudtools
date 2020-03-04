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
  cout << "Convert a point cloud and trajectory file into a ray cloud" << endl;
  cout << "usage:" << endl;
  cout << "rayconvert pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << endl;
  cout << "                                            trajectoryfile is a text file in time,x,y,z,qw,qx,qy,qz format" << endl;
  cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();
  
  string pointCloudFile = argv[1];
  string trajectoryFile = argv[2];
  Cloud cloud;
  cloud.load(pointCloudFile, trajectoryFile);

  string saveFile = pointCloudFile.substr(0, pointCloudFile.size()-4);
  if (pointCloudFile.substr(pointCloudFile.size()-4)==".ply")
    saveFile += "_raycloud";
  cloud.save(saveFile + ".ply");

  return true;
}
