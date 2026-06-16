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
#include <limits>

#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raydiffer.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Difference between two ray clouds output a coloured cloud, cloud1 differences in scarlet, cloud2 differences in green/cyan, and similarity printed to screen." << std::endl;
  std::cout << "usage (--view / -v to view results):" << std::endl;
  std::cout << "raydiff cloud1.ply cloud2.ply" << std::endl;
  std::cout << "                              --distance 0 - optional threshold in m for colouring differences. Default auto-detects distribution shoulder" << std::endl;
  std::cout << "                              --individual - output as two clouds, to show differences on each." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayDiff(int argc, char *argv[])
{
  ray::FileArgument cloud1_name, cloud2_name;
  ray::OptionalFlagArgument visualise("view", 'v'), individual("individual", 'i');
  ray::DoubleArgument distance_threshold(0.0, 1000.0, 0.0);
  ray::OptionalKeyValueArgument distance_option("distance", 'd', &distance_threshold);
  if (!ray::parseCommandLine(argc, argv, { &cloud1_name, &cloud2_name }, { &distance_option, &visualise, &individual }))
  {
    usage();
  }

  // Next I need to convert raydiff into a maplib function so I can call it efficiently
  // 1. the input needn't be a ray cloud, just a list of points...
  // 2. the output doesn't have to be a ray cloud

  // low memory raydiff requires offsets to support float vectors, and careful avoidance of duplicating large arrays
  // step 1: get points only
  std::vector<Eigen::Vector3f> points1, points2;
  std::vector<Eigen::Vector3f> *points = &points1;
  Eigen::Vector3d start_pos(0,0,0);
  bool has_start_pos = false;
  auto getPoints = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &, std::vector<ray::RGBA> &colours) 
  {
    if (!has_start_pos && !ends.empty())
    {
      start_pos = ends[0];
      has_start_pos = true;
    }
    for (size_t i = 0; i<ends.size(); i++)
    {
      if (colours[i].alpha == 0)
        points->push_back(Eigen::Vector3f(std::numeric_limits<float>::max(), 0,0)); // flag as unbounded
      else
        points->push_back((ends[i] - start_pos).cast<float>());
    }
  };

  if (!ray::Cloud::read(cloud1_name.name(), getPoints))
    return false;  
  points = &points2;
  if (!ray::Cloud::read(cloud2_name.name(), getPoints))
    return false;  
  
  
  // Primary goal - get an overall similarity percentage and return it as the return value, so it can be used in bash scripts
  // this needs to include colour and location, but be insensitive to changes in density. i.e. it needs to work usefully on a repeat scan.
  std::vector<float> dists_to_cloud2, dists_to_cloud1;
  ray::getDistancesBetweenPoints(points1, points2, dists_to_cloud1, dists_to_cloud2);
  points1.clear(); points1.shrink_to_fit();
  points2.clear(); points2.shrink_to_fit();

  double dist_threshold = distance_threshold.value();
  double similarity = ray::printDistanceStatistics(dists_to_cloud1, dists_to_cloud2, dist_threshold);
  
  if (!ray::writeDifferencesToRayClouds(cloud1_name.nameStub(), cloud2_name.nameStub(), dists_to_cloud1, dists_to_cloud2, dist_threshold, individual.isSet(), visualise.isSet()))
  {
    return 1;
  }

  return (int)similarity;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDiff, argc, argv);
}