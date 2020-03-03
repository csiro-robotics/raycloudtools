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
  cout << "Decimate a ray cloud spatially or temporally" << endl;
  cout << "usage:" << endl;
  cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << endl;
  cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << endl;
  cout << "raydecimate raycloud 0.1 seconds - reduces to one ray every 0.1 seconds" << endl;
  exit(error);
}

#if 0 // could be useful for decimate?
struct Vector3iLess
{
  bool operator()(const Vector3i &a, const Vector3i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    return a[2] < b[2];
  }
};

vector<int> voxelRandomSubsample(const LaserData &laser, const vector<int> &subset, double voxelWidth)
{
  vector<int> indices;
  set<Vector3i, Vector3iLess> testSet;
  for (unsigned int j = 0; j<subset.size(); j++)
  {
    int i = subset[j];
    Vector3i places = Vector3i(floor(laser.positions[i][0] / voxelWidth), floor(laser.positions[i][1] / voxelWidth), floor(laser.positions[i][2] / voxelWidth));
    if (testSet.find(places) == testSet.end())
    {
      testSet.insert(places);
      indices.push_back(i);
    }
  }
  return indices;
}
#endif

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4)
    usage();
  
  string file = argv[1];
  Cloud cloud;
  cloud.load(file);
  string type = argv[3];
  if (type=="cm")
  {
    cloud.decimate(stod(argv[2]));
  }
  else if (type == "rays")
  {
    int decimation = stoi(argv[2]);
    vector<Vector3d> starts, ends, times, intensities;
    for (int i = 0; i<cloud.ends.size(); i+=decimation)
    {
      starts.push_back(cloud.starts[i]);
      ends.push_back(cloud.ends[i]);
      times.push_back(cloud.times[i]);
      intensities.push_back(cloud.intensities[i]);
    }
    cloud.starts = starts; cloud.ends = ends; cloud.times = times; cloud.intensities = intensities;
  }
  else if (type == "seconds" || type == "s")
  {
    double delta = stod(argv[2]);
    vector<Vector3d> starts, ends, times, intensities;
    double lastTime = cloud.times[0];
    for (int i = 0; i<cloud.ends.size(); i++)
    {
      if (cloud.times[i] >= lastTime + delta)
      {
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        intensities.push_back(cloud.intensities[i]);
        while (cloud.times[i] >= lastTime + delta) // in case delta is tiny
          lastTime += delta;
      }
    }
    cloud.starts = starts; cloud.ends = ends; cloud.times = times; cloud.intensities = intensities;
  }
  else 
    usage(false);
  cloud.save(file);
  return true;
}
