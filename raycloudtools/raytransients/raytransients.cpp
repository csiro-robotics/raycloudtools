// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymerger.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raylib/raythreads.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <thread>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Splits a raycloud into the transient rays and the fixed part" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytransients min raycloud 20 rays - splits out positive transients (objects that have since moved)." << std::endl;
  std::cout << "                                     20 is number of pass through rays to classify as transient." << std::endl;
  std::cout << "              max    - finds negative transients, such as a hallway exposed when a door opens." << std::endl;
  std::cout << "              oldest - keeps the oldest geometry when there is a difference over time." << std::endl;
  std::cout << "              newest - uses the newest geometry when there is a difference over time." << std::endl;
  std::cout << " --colour     - also colours the clouds, to help tweak numRays. blue: opacity, green: pass throughs." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayTransients(int argc, char *argv[])
{
  ray::KeyChoice merge_type({ "min", "max", "oldest", "newest" });
  ray::FileArgument cloud_file;
  ray::DoubleArgument num_rays(0.1, 100.0);
  ray::TextArgument text("rays");
  ray::OptionalFlagArgument colour("colour", 'c');
  if (!ray::parseCommandLine(argc, argv, { &merge_type, &cloud_file, &num_rays, &text }, { &colour }))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();

  ray::Threads::init();
  ray::MergerConfig config;
  // Note: we actually get better multi-threaded performace with smaller voxels
  config.voxel_size = 0.0;
  config.num_rays_filter_threshold = num_rays.value();
  config.merge_type = ray::MergeType::Mininum;
  config.colour_cloud = colour.isSet();

  if (merge_type.selectedKey() == "oldest")
  {
    config.merge_type = ray::MergeType::Oldest;
  }
  if (merge_type.selectedKey() == "newest")
  {
    config.merge_type = ray::MergeType::Newest;
  }
  if (merge_type.selectedKey() == "min")
  {
    config.merge_type = ray::MergeType::Mininum;
  }
  if (merge_type.selectedKey() == "max")
  {
    config.merge_type = ray::MergeType::Maximum;
  }

  ray::Merger filter(config);
  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);

  filter.filter(cloud, &progress);

  progress_thread.requestQuit();
  progress_thread.join();

  const ray::Cloud &transient = filter.differenceCloud();
  const ray::Cloud &fixed = filter.fixedCloud();

  transient.save(cloud_file.nameStub() + "_transient.ply");
  fixed.save(cloud_file.nameStub() + "_fixed.ply");
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayTransients, argc, argv);
}