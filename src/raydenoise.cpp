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
  cout << "Remove noise from ray clouds. In particular edge noise and isolated point noise." << endl;
  cout << "usage:" << endl;
  cout << "raydenoise raycloud 3 cm   - removes rays that contact more than 3 cm from any other," << endl;
  cout << "                             or are more than 3cm different in range from temporally adjacent rays." << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4)
    usage();
  if (argv[3] != "cm")
    usage();
  string file = argv[1];
  Cloud cloud;
  cloud.load(file);
  double distance = stod(argv[2]);

  // Firstly look at adjacent rays by range. We don't want to throw away large changes,
  // instead, the intermediate of 3 adjacent ranges that is too far from both ends...
  Cloud newCloud;
  for (int i = 1; i<(int)cloud.rayStarts.size()-1; i++)
  {
    double range0 = (cloud.rayEnds[i-1]-cloud.rayStarts[i-1]).norm();
    double range1 = (cloud.rayEnds[i]-cloud.rayStarts[i]).norm();
    double range2 = (cloud.rayEnds[i+1]-cloud.rayStarts[i+1]).norm();
    double minDist = min(abs(range1-range0), abs(range2 - range0));
    if (minDist > distance)
    {
      cloud.starts.erase(cloud.starts.begin()+i);
      cloud.ends.erase(cloud.ends.begin()+i);
      cloud.times.erase(cloud.times.begin()+i);
      cloud.intensities.erase(cloud.intensities.begin()+i);
    }

    // Next, look at spatially isolated points:
    // TODO: requires NABO...
  }

  return true;
}
