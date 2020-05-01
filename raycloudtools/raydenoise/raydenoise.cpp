// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace ray;

void usage(int exit_code = 0)
{
  cout << "Remove noise from ray clouds. In particular edge noise and isolated point noise." << endl;
  cout << "usage:" << endl;
  cout << "raydenoise raycloud 3 cm   - removes rays that contact more than 3 cm from any other," << endl;
  cout << "           --range 4 cm    - option to also remove mixed-signal noise that occurs at a range gap." << endl;
  exit(exit_code);
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
    double range_distance = 0.01 * stod(argv[5]);
    if (range_distance < 0.01)
      usage();
    
    Cloud new_cloud;
    new_cloud.starts.reserve(cloud.starts.size());
    new_cloud.ends.reserve(cloud.ends.size());
    new_cloud.times.reserve(cloud.times.size());
    new_cloud.colours.reserve(cloud.colours.size());
    // Firstly look at adjacent rays by range. We don't want to throw away large changes,
    // instead, the intermediate of 3 adjacent ranges that is too far from both ends...
    for (int i = 1; i<(int)cloud.starts.size()-1; i++)
    {
      double range0 = (cloud.ends[i-1]-cloud.starts[i-1]).norm();
      double range1 = (cloud.ends[i]-cloud.starts[i]).norm();
      double range2 = (cloud.ends[i+1]-cloud.starts[i+1]).norm();
      double min_dist = min(abs(range0 - range2), min(abs(range1-range0), abs(range2 - range1)));
      if (!cloud.rayBounded(i) || min_dist < range_distance)
      {
        new_cloud.starts.push_back(cloud.starts[i]);
        new_cloud.ends.push_back(cloud.ends[i]);
        new_cloud.times.push_back(cloud.times[i]);
        new_cloud.colours.push_back(cloud.colours[i]);
      }
    }
    cout << cloud.starts.size()-new_cloud.starts.size() << " rays removed with range gaps > " << range_distance*100.0 << " cm." << endl;
    cloud = new_cloud;
  }

  // Next, look at spatially isolated points:
  {
    MatrixXi indices;
    MatrixXd dists2;

    // simplest scheme... find 3 nearest neighbours and do cross product
    Nabo::NNSearchD *nns;
    Nabo::Parameters params("bucketSize", 8);
    vector<Vector3d> &points = cloud.ends;
    MatrixXd points_p(3, points.size());
    for (unsigned int i = 0; i<points.size(); i++)
      points_p.col(i) = points[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);//, 0, params);

    // Run the search
    const int search_size = 10;
    indices.resize(search_size, points.size());
    dists2.resize(search_size, points.size());
    nns->knn(points_p, indices, dists2, search_size, 0.01, 0, 1.0);
    delete nns;

    Cloud new_cloud;
    new_cloud.starts.reserve(cloud.starts.size());
    new_cloud.ends.reserve(cloud.ends.size());
    new_cloud.times.reserve(cloud.times.size());
    new_cloud.colours.reserve(cloud.colours.size());
    for (int i = 0; i<(int)points.size(); i++)
    {
      if (!cloud.rayBounded(i) || (dists2(0,i) < 1e10 && dists2(0,i) < sqr(distance)))
      {
        new_cloud.starts.push_back(cloud.starts[i]);
        new_cloud.ends.push_back(cloud.ends[i]);
        new_cloud.times.push_back(cloud.times[i]);
        new_cloud.colours.push_back(cloud.colours[i]);
      }
    }
    cout << cloud.starts.size()-new_cloud.starts.size() << " rays removed with ends further than " << distance*100.0 << " cm from any other." << endl;
    cloud = new_cloud;
  }

  cloud.save(file.substr(0,file.length()-4) + "_denoised.ply");

  return true;
}
