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
#include <limits>

#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/rayply.h"
#include "raylib/raycloudwriter.h"
#include "raylib/raydiffer.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Difference between two ray clouds output a coloured cloud, cloud1 differences in red, cloud2 differences in green, and similarity printed to screen." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydiff cloud1.ply cloud2.ply" << std::endl;
  std::cout << "                              --distance 0 - optional threshold in m for colouring differences. Default auto-detects distribution shoulder" << std::endl;
  std::cout << "                              --individual - output as two clouds, to show differences on each." << std::endl;
  std::cout << "                              --visualise  - open in the default visualisation tool" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayDiff(int argc, char *argv[])
{
  ray::FileArgument cloud1_name, cloud2_name;
  ray::OptionalFlagArgument visualise("visualise", 'v'), individual("individual", 'i');
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

  std::cout << "saving out coloured red for matches beyond shoulder difference:" << std::endl;
  
  // now render visuals
  ray::CloudWriter writer;
  if (!writer.begin(cloud1_name.nameStub() + "_diff.ply"))
    return false;

  // By maintaining these buffers below, we avoid almost all memory fragmentation
  ray::Cloud chunk;
  std::map<Eigen::Vector3i, Eigen::Vector2i, ray::Vector3iLess> voxel_map;
  std::vector<Eigen::Vector3i> samples;

  int j = 0;
  Eigen::Vector3d diff_col(255,0,0);
  std::vector<float> *dists = &dists_to_cloud1;
  const float eps = 1e-8f; // in case there is inaccuracy in the KNN distance estimation for co-located point pairs
  auto colour = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                    std::vector<double> &times, std::vector<ray::RGBA> &colours) 
  {
    // firstly we store a count per cell
    chunk.clear();
    for (size_t i = 0; i<ends.size(); i++)
    {
      if (colours[i].alpha > 0)
      {
        if ((*dists)[j] > eps)
        {
          chunk.addRay(starts[i], ends[i], times[i], colours[i]);
          if ((*dists)[j] > dist_threshold)
          {
            double proximity = dist_threshold / (*dists)[j];
            Eigen::Vector3d col(colours[i].red, colours[i].green, colours[i].blue);
            Eigen::Vector3d new_col = diff_col + (col - diff_col) * proximity;
            chunk.colours.back() = ray::RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), colours[i].alpha);
          }
        }
        j++;
      }
      else
        chunk.addRay(starts[i], ends[i], times[i], colours[i]);
    }
    writer.writeChunk(chunk);
  };

  if (!ray::Cloud::read(cloud1_name.name(), colour))
    return false;
  j = 0;
  dists_to_cloud1.clear();
  dists_to_cloud1.shrink_to_fit();
  Eigen::Vector3d diff2_col(0,255,0);

  dists = &dists_to_cloud2;
  if (!individual.isSet())
  {
    auto colour2 = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours) 
    {
      // firstly we store a count per cell
      chunk.clear();
      for (size_t i = 0; i<ends.size(); i++)
      {
        if (colours[i].alpha > 0)
        {
          if ((*dists)[j] > eps)
          {
            chunk.addRay(starts[i], ends[i], times[i], colours[i]);
            if ((*dists)[j] > dist_threshold)
            {
              double proximity = dist_threshold / (*dists)[j];
              Eigen::Vector3d col(colours[i].red, colours[i].green, colours[i].blue);
              Eigen::Vector3d new_col = diff2_col + (col - diff2_col) * proximity;
              chunk.colours[i] = ray::RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), colours[i].alpha);
            }
          }
          j++;
        }
        else
          chunk.addRay(starts[i], ends[i], times[i], colours[i]);
      }
      writer.writeChunk(chunk);
    };

    if (!ray::Cloud::read(cloud2_name.name(), colour2))
      return false;
    writer.end();

    if (visualise.isSet())
    {
      std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_name.nameStub() + "_diff.ply";
      system(command.c_str());  
    }
  }
  else
  {
    writer.end();
    diff_col = diff2_col;
    if (!writer.begin(cloud2_name.nameStub() + "_diff.ply"))
      return false;

    if (!ray::Cloud::read(cloud2_name.name(), colour))
      return false;
    writer.end();
    if (visualise.isSet())
    {
      std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_name.nameStub() + "_diff.ply " + cloud2_name.nameStub() + "_diff.ply";
      system(command.c_str());  
    }
  }

  return (int)similarity;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDiff, argc, argv);
}