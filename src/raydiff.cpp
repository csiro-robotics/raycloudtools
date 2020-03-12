// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayutils.h"
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
  cout << "Reports the difference between two rayclouds." << endl;
  cout << "usage:" << endl;
  cout << "raydiff raycloud1 raycloud2 - name of the two ray clouds" << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc != 3)
    usage();
  
  string file1 = argv[1];
  string file2 = argv[2];
  Cloud cloud1, cloud2;
  cloud1.load(file1);
  cloud2.load(file2);
#if 0
  // Firstly get the transform of best fit between the clouds
  double timeShift = 0.0;
  Pose poseChange = getTransformOfBestFit(cloud1, cloud2, timeShift);
  cout << "Difference in pose: " << poseChange.position.transpose() << ", " << poseChange.rotation << endl;
  // and apply the transform to cloud 2
  cloud2.transform(poseChange);

  // Next, estimate the overlap volume as a grid...
  Grid grid;
  grid.estimateOverlap(cloud1, cloud2);
  int numOccupied = 0;
  int numSimilar = 0;
  for (auto &cell: cells)
  {
    if (cell.ends.size()>5)
      numOccupied++;
    if (cell.overlaps)
      numSimilar++;
  }
  double overlap = 100.0*(double)numSimilar/numOccupied;
  cout << "overlap percentage: " << overlap << endl;
  // and crop out the non-overlapping parts
  for (auto &cell: cells)
  {
    if (!cell.overlaps)
    {
      cloud1.crop(cell.bounds);
      cloud2.crop(cell.bounds);
    }
  }

  // Now, find the relative error between the two clouds.
  double error = mean(nearestNeighbourDistances(cloud1, cloud2));
  cout << "mean difference in overlapping region: " << error << " metres" << endl;

  // Lastly, summarise the differences
  int numErrors = 0;
  string errors;
  if (poseChange.position.norm() > 1.0)
  {
    numErrors++;
    string += "translation > 1m, ";
  }
  if (posehange.rotation.norm() > 0.1)
  {
    numErrors++;
    string += "rotation quat > 0.1, ";
  }
  if (abs(timeShift) > 1.0)
  {
    numErrors++;
    string += "time > 1s, ";
  }
  if (overlap < 90)
  {
    numErrors++;
    string += "overlap < 90\%, ";
  }
  if (error > 0.1)
  {
    numErrors++;
    string += "difference inside overlap > 0.1m, "
  }
  cout << endl;
  cout << numErrors << " key differences." << endl;
  cout << "(" << substr(string, 0, string.length()-2) << ")." << endl;

#endif
  return true;
}
