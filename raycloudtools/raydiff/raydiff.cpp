// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <chrono>
#include <nabo/nabo.h>

#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/rayply.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Difference between two ray clouds. Optional visualisation." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydiff cloud1.ply cloud2.ply" << std::endl;
  std::cout << "                              --visualise     - open in the default visualisation tool" << std::endl;
  // clang-format on
  exit(exit_code);
}

void calcNearestNeighbourDistances(const ray::Cloud &cloud1, const ray::Cloud &cloud2, std::vector<double> &distances)
{
  Nabo::NNSearchD *nns;
 // Nabo::Parameters params("bucketSize", 8);
  int num_bounded = 0, num_bounded2 = 0;
  for (unsigned int i = 0; i < cloud1.ends.size(); i++)
    if (cloud1.rayBounded(i))
      num_bounded++;
  for (unsigned int i = 0; i < cloud2.ends.size(); i++)
    if (cloud2.rayBounded(i))
      num_bounded2++;
  Eigen::MatrixXd points_p(3, num_bounded);
  int j = 0;
  for (unsigned int i = 0; i < cloud1.ends.size(); i++) 
  {
    if (cloud1.rayBounded(i))
      points_p.col(j++) = cloud1.ends[i];
  }
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  Eigen::MatrixXd points_q(3, num_bounded2);
  j = 0;
  for (unsigned int i = 0; i < cloud2.ends.size(); i++) 
  {
    if (cloud2.rayBounded(i))
      points_q.col(j++) = cloud2.ends[i];
  }
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  // Run the search
  const int search_size = 1;
  indices.resize(search_size, num_bounded2);
  dists2.resize(search_size, num_bounded2);
  nns->knn(points_q, indices, dists2, search_size, ray::kNearestNeighbourEpsilon, Nabo::NNSearchD::SearchOptionFlags::SORT_RESULTS | Nabo::NNSearchD::SearchOptionFlags::ALLOW_SELF_MATCH);
  distances.resize(num_bounded2);
  for (int i = 0; i<num_bounded2; i++)
  {
    distances[i] = std::sqrt(dists2(0,i));
  }
  delete nns;
}

int rayDiff(int argc, char *argv[])
{
  ray::FileArgument cloud1_name, cloud2_name;
  ray::OptionalFlagArgument visualise("--visualise", 'v');
  if (!ray::parseCommandLine(argc, argv, { &cloud1_name, &cloud2_name }, { &visualise }))
  {
    usage();
  }

  ray::Cloud cloud1, cloud2;
  if (!cloud1.load(cloud1_name.name()) || !cloud2.load(cloud2_name.name()))
  {
    usage();
  }

  // Primary goal - get an overall similarity percentage and return it as the return value, so it can be used in bash scripts
  // this needs to include colour and location, but be insensitive to changes in density. i.e. it needs to work usefully on a repeat scan.
  std::vector<double> dists_to_cloud2, dists_to_cloud1;
  calcNearestNeighbourDistances(cloud1, cloud2, dists_to_cloud2);
  calcNearestNeighbourDistances(cloud2, cloud1, dists_to_cloud1);

  // 1. max distance
  double max_dist = 0.0;
  // 2. percentage of cloud within tolerance

  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud2[i]);
  for (int i = 0; i<(int)dists_to_cloud1.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud1[i]);

  // 3. median distance
  dists_to_cloud2.insert(dists_to_cloud2.end(), dists_to_cloud1.begin(), dists_to_cloud1.end());
  double median_dist = ray::median(dists_to_cloud2);

  std::sort(dists_to_cloud2.begin(), dists_to_cloud2.end());
  double cum_dist = 0.0;
  double cum_area = 0.0;
  double last_dist = 0.0;
  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
  {
    double dist = dists_to_cloud2[i];
    cum_dist += dist;
    double density = cum_dist / dist;

    cum_area += density * (dist - last_dist);
    last_dist = dist; 
  }
  double total_area = cum_area;
  
  cum_dist = 0.0;
  cum_area = 0.0;
  last_dist = 0.0;
  double shoulder_error = 0.0;
  double shoulder_uniformity = 0.0;
  double last_dif = 0.0;
  int shoulder_i = 0;
  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
  {
    double dist = dists_to_cloud2[i];
    double count = (double)(i+1);
    double density = count / dist;
    cum_area += density * (dist - last_dist);
    double dif = (cum_area - count) - (total_area - cum_area);
    if (dif >= 0.0)
    {
      // linearly interpolate between entries
      double blend = -last_dif / (dif - last_dif);
      shoulder_error = last_dist + (dist - last_dist)*blend;
      shoulder_i = i;
      
      shoulder_uniformity = ((count - 1.0) + blend) / total_area;
      break;
    }
    last_dif = dif;
    last_dist = dist; 
  }
  std::cout << "median difference: " << median_dist << " m" << std::endl;
  std::cout << "max difference: " << max_dist << " m" << std::endl;
  std::cout << "shoulder difference: " << shoulder_error << " m" << std::endl;
  std::cout << "inside shoulder: " << 100.0*(double)shoulder_i / (double)dists_to_cloud2.size() << "%" << std::endl;
  std::cout << "shoulder uniformity: " << 100.0*shoulder_uniformity << "%" << std::endl;

  if (visualise.isSet())
  {
    std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_name.nameStub() + ".ply";
    return system(command.c_str());  
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDiff, argc, argv);
}