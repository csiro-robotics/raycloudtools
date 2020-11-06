// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits>
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raylaz.h"
#include "raylib/rayply.h"
#include "raylib/raytrajectory.h"

void usage(int exit_code = 1)
{
  std::cout << "Export a ray cloud into a point cloud amd trajectory file" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayexport raycloudfile.ply pointcloud.ply trajectoryfile.ply - output in specified formats"
       << std::endl;
  std::cout << "                           pointcloud.laz trajectoryfile.txt" << std::endl;
  std::cout << "                           --frequency 100 - maximum trajectory frequency (Hz). Default is all points" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument raycloud_file, pointcloud_file, trajectory_file;
  ray::DoubleArgument frequency(0.0001, 10000);
  ray::OptionalKeyValueArgument frequency_option("frequency", 'f', &frequency);
  if (!ray::parseCommandLine(argc, argv, {&raycloud_file, &pointcloud_file, &trajectory_file}, {&frequency_option}))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(raycloud_file.name()))
    usage();

  if (pointcloud_file.nameExt() == "laz")
    writeLas(pointcloud_file.name(), cloud.ends, cloud.times, cloud.colours);
  else if (pointcloud_file.nameExt() == "ply")
    writePlyPointCloud(pointcloud_file.name(), cloud.ends, cloud.times, cloud.colours);
  else
    usage();

  // In order to temporally decimate the trajectory, we have to make sure it is sorted first.
  struct TimeIndex
  {
    double time;
    size_t index;
  };
  std::vector<TimeIndex> time_indices(cloud.times.size());
  bool sorted = true;
  for (size_t i = 0; i<cloud.times.size(); i++)
  {
    time_indices[i].time = cloud.times[i];
    time_indices[i].index = i;
    if (i > 0 && cloud.times[i] < cloud.times[i-1])
      sorted = false;
  }
  if (!sorted)
    std::sort(time_indices.begin(), time_indices.end(), [](const TimeIndex &a, const TimeIndex &b){ return a.time < b.time; });

  // now temporally decimate:
  std::vector<Eigen::Vector3d> starts;
  std::vector<double> times;
  std::vector<ray::RGBA> colours;
  double time_step = 1.0 / frequency.value();
  if (frequency_option.isSet())
  {
    double last_time = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i<time_indices.size(); i++)
    {
      if (time_indices[i].time >= last_time + time_step)
      {
        size_t id = time_indices[i].index;
        starts.push_back(cloud.starts[id]);
        times.push_back(cloud.times[id]);
        colours.push_back(cloud.colours[id]);
        last_time = time_indices[i].time;
      }
    }
  }
  else // apply the sorted trajectory, even if we're not using a frequency
  {
    for (size_t i = 0; i<time_indices.size(); i++)
    {
      size_t id = time_indices[i].index;
      starts.push_back(cloud.starts[id]);
      times.push_back(cloud.times[id]);
      colours.push_back(cloud.colours[id]);
    }
  }

  // and output the result, depending on the file format
  if (trajectory_file.nameExt() == "txt")
  {
    ray::Trajectory trajectory;
    trajectory.times() = std::move(times);
    trajectory.points() = std::move(starts);
    trajectory.save(trajectory_file.name());
  }
  else if (trajectory_file.nameExt() == "ply")
    writePlyPointCloud(trajectory_file.name(), starts, times, colours);
  else
    usage();
}      
