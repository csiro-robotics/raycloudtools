#include "raycloud.h"
#include <nabo/nabo.h>

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
  cout << "           --range 4 cm    - option to also remove mixed-signal noise that occurs at a range gap." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 7)
    usage();
  if (string(argv[3]) != "cm")
    usage();

  string file = argv[1];
  Cloud cloud;
  cloud.load(file);
  double distance = 0.01 * stod(argv[2]);
  if (distance < 0.01)
    usage();
  
  // Look at gaps in range... this signals a mixed-signal where the lidar has contacted two surfaces.
  if (argc == 7)
  {
    if ((string(argv[4]) != "--range" && string(argv[4]) != "-r") || string(argv[6]) != "cm")
      usage();
    double rangeDistance = 0.01 * stod(argv[5]);
    if (rangeDistance < 0.01)
      usage();
    
    vector<Vector3d> starts, ends;
    vector<double> times, intensities;
    starts.reserve(cloud.starts.size());
    ends.reserve(cloud.ends.size());
    times.reserve(cloud.times.size());
    intensities.reserve(cloud.intensities.size());
    // Firstly look at adjacent rays by range. We don't want to throw away large changes,
    // instead, the intermediate of 3 adjacent ranges that is too far from both ends...
    for (int i = 1; i<(int)cloud.starts.size()-1; i++)
    {
      double range0 = (cloud.ends[i-1]-cloud.starts[i-1]).norm();
      double range1 = (cloud.ends[i]-cloud.starts[i]).norm();
      double range2 = (cloud.ends[i+1]-cloud.starts[i+1]).norm();
      double minDist = min(abs(range0 - range2), min(abs(range1-range0), abs(range2 - range1)));
      if (minDist < rangeDistance)
      {
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        if (cloud.intensities.size() > 0)
          intensities.push_back(cloud.intensities[i]);
      }
    }
    cout << cloud.starts.size()-starts.size() << " rays removed with range gaps > " << rangeDistance*100.0 << " cm." << endl;
    cloud.starts = starts;
    cloud.ends = ends;
    cloud.times = times;
    cloud.intensities = intensities;
  }

  // Next, look at spatially isolated points:
  {
    MatrixXi indices;
    MatrixXd dists2;

    // simplest scheme... find 3 nearest neighbours and do cross product
    Nabo::NNSearchD *nns;
    Nabo::Parameters params("bucketSize", 8);
    vector<Vector3d> &points = cloud.ends;
    MatrixXd pointsP(3, points.size());
    for (unsigned int i = 0; i<points.size(); i++)
      pointsP.col(i) = points[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);//, 0, params);

    // Run the search
    const int searchSize = 10;
    indices.resize(searchSize, points.size());
    dists2.resize(searchSize, points.size());
    nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0);
    delete nns;

    vector<char> ignoreList(points.size());
    memset(&ignoreList[0], 0, sizeof(bool) * ignoreList.size());
    int ignored = 0;

    vector<Vector3d> starts, ends;
    vector<double> times, intensities;
    starts.reserve(cloud.starts.size());
    ends.reserve(cloud.ends.size());
    times.reserve(cloud.times.size());
    intensities.reserve(cloud.intensities.size());
    for (int i = 0; i<points.size(); i++)
    {
      if (dists2(0,i) < 1e10 && dists2(0,i) < sqr(distance))
      {
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        if (cloud.intensities.size() > 0)
          intensities.push_back(cloud.intensities[i]);
      }
    }
    cout << cloud.starts.size()-starts.size() << " rays removed with ends further than " << distance*100.0 << " cm from any other." << endl;
    cloud.starts = starts;
    cloud.ends = ends;
    cloud.times = times;
    cloud.intensities = intensities;
  }

  cloud.save(file.substr(0,file.length()-4) + "_denoised.ply");

  return true;
}
