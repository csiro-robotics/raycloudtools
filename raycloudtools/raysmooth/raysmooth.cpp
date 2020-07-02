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

void usage(int exit_code = 0)
{
  std::cout << "Smooth a ray cloud. Nearby off-surface points are moved onto the nearest surface." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysmooth raycloud" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 2)
    usage();

  std::string file = argv[1];
  ray::Cloud cloud;
  cloud.load(file);

  // Method:
  // 1. generate normals and neighbour indices
  // 2. pull point along normal direction so as to match neighbours, weighted by normal similarity

  const int num_neighbours = 16;
  std::vector<Eigen::Vector3d> normals;
  Eigen::MatrixXi neighbour_indices;
  cloud.getSurfels(num_neighbours, NULL, &normals, NULL, NULL, &neighbour_indices);

  std::vector<Eigen::Vector3d> centroids(cloud.ends.size());
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    double total_weight = 0.2; // more averaging if it uses less of the central position, but 0 risks a divide by 0
    Eigen::Vector3d weighted_sum = cloud.ends[i]*total_weight;
    for (int j = 0; j < num_neighbours && neighbour_indices(j, i) > -1; j++) 
    {
      int k = neighbour_indices(j, i);
      double weight = std::max(0.0, 1.0 - (normals[k] - normals[i]).squaredNorm());
      weighted_sum += cloud.ends[k] * weight;
      total_weight += weight;
    }
    centroids[i] = weighted_sum / total_weight;
  }
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    cloud.ends[i] += normals[i] * (centroids[i]-cloud.ends[i]).dot(normals[i]);
  }

  cloud.save(file.substr(0, file.length() - 4) + "_smooth.ply");

  return true;
}
