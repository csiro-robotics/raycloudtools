// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout << "Decimate a ray cloud spatially or temporally" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << std::endl;
  std::cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::IntArgument num_rays(1,100);
  ray::DoubleArgument vox_width(0.01, 100.0);
  ray::ValueKeyChoice quantity({&vox_width, &num_rays}, {"cm", "rays"});
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &quantity}))
    usage();

  std::ofstream ofs;
  if (!ray::writePlyChunkStart(cloud_file.nameStub() + "_decimated.ply", ofs))
    usage();

  // By maintaining these buffers below, we avoid almost all memory fragmentation  
  std::vector<Eigen::Matrix<float, 9, 1>> buffer;
  std::vector<Eigen::Vector3d> start_points;
  std::vector<Eigen::Vector3d> end_points;
  std::vector<double> time_points;
  std::vector<ray::RGBA> colour_values;
  std::vector<int64_t> subsample;

  // voxel set is global, however its size is proportional to the decimated cloud size,
  // so we expect it to fit within RAM limits
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;  

  auto decimate_cm = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    double width = 0.01 * vox_width.value();

    subsample.clear();
    voxelSubsample(ends, width, subsample, &voxel_set);
    start_points.resize(subsample.size());
    end_points.resize(subsample.size());
    time_points.resize(subsample.size());
    colour_values.resize(subsample.size());
    for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
    {
      int64_t id = subsample[i];
      start_points[i] = starts[id];
      end_points[i] = ends[id];
      colour_values[i] = colours[id];
      time_points[i] = times[id];
    }
    ray::writePlyChunk(ofs, buffer, start_points, end_points, time_points, colour_values);
  };
  auto decimate_rays = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    int decimation = num_rays.value();
    size_t count = 1 + ends.size() / (size_t)decimation;
    start_points.resize(count);
    end_points.resize(count);
    time_points.resize(count);
    colour_values.resize(count);
    for (size_t i = 0, c = 0; i < ends.size(); i += decimation, c++)
    {
      start_points[c] = starts[i];
      end_points[c] = ends[i];
      time_points[c] = times[i];
      colour_values[c] = colours[i];
    }
    ray::writePlyChunk(ofs, buffer, start_points, end_points, time_points, colour_values);
  };

  bool success;
  if (quantity.selectedKey() == "cm")
    success = ray::readPly(cloud_file.name(), false, decimate_cm);
  else
    success = ray::readPly(cloud_file.name(), false, decimate_rays);
  if (!success)
    usage();
  ray::writePlyChunkEnd(ofs);

  return 0;
}
