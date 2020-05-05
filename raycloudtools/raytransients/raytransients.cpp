// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raydebugdraw.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayprogress.h"
#include "raylib/raytransientfilter.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <thread>

using namespace std;
using namespace Eigen;
using namespace ray;

#define NEW_FILTER 1

void usage(int exit_code = 0)
{
  cout << "Splits a raycloud into the transient rays and the fixed part" << endl;
  cout << "usage:" << endl;
  cout << "raytransients min raycloud 20 rays - splits out positive transients (objects that have since moved)."
       << endl;
  cout << "                                     20 is number of pass through rays to classify as transient." << endl;
  cout << "              max    - finds negative transients, such as a hallway exposed when a door opens." << endl;
  cout << "              oldest - keeps the oldest geometry when there is a difference over time." << endl;
  cout << "              newest - uses the newest geometry when there is a difference over time." << endl;
  cout << " --colour     - also colours the clouds, to help tweak numRays. red: opacity, green: pass throughs, blue: "
          "planarity."
       << endl;
  exit(exit_code);
}

void runProrgess(const Progress &progress, std::atomic_bool &quit)
{
  Progress last;
  Progress current;
  progress.read(&last);

  const auto show_progress = [] (Progress &p, bool finalise) //
  {
    if (finalise)
    {
      // Finalise progress == target.
      if (p.target())
      {
        p.setProgress(p.target());
      }
    }

    if (p.phase().length() || p.target() || p.progress())
    {
      std::cout << "\r                                    \r";
      std::cout << p.phase() << ' ' << p.progress();
      if (size_t target = p.target())
      {
        std::cout << " / " << target;
      }

      if (finalise)
      {
        std::cout << std::endl;
      }
      else
      {
        std::cout << std::flush;
      }
    }
  };

  while (!quit)
  {
    progress.read(&current);
    if (current.phase() != last.phase())
    {
      // Ensure we finalise the display.
      show_progress(last, true);
    }

    if ( current.progress() != last.progress() || current.target() != last.target())
    {
      show_progress(current, false);
      current.read(&last);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }

  // Past update.
  progress.read(&current);
  // Do not finalise in case we didn't go full progress.
  show_progress(last, false);
  std::cout << std::endl;
}

int main(int argc, char *argv[])
{
  DebugDraw::init(argc, argv, "raytransients");
  if (argc != 5 && argc != 6)
    usage();

  bool colour = false;
  if (argc == 6)
  {
    if (string(argv[5]) != "--colour" && string(argv[5]) != "-c")
      usage();
    colour = true;
  }
  double num_rays = stod(argv[3]);
  string merge_type = argv[1];
  if (merge_type != "min" && merge_type != "max" && merge_type != "oldest" && merge_type != "newest")
    usage();
  string file = argv[2];
  Cloud cloud;
  cloud.load(file);

#if NEW_FILTER
  TransientFilterConfig config;
  config.merge_type = TransientFilterType::Maximum;
  config.num_rays_filter_threshold = num_rays;
  config.colour_cloud = colour;
  config.voxel_size = 0.25;

  if (merge_type == "oldest")
  {
    config.merge_type = TransientFilterType::Oldest;
  }
  if (merge_type == "newest")
  {
    config.merge_type = TransientFilterType::Newest;
  }
  if (merge_type == "min")
  {
    config.merge_type = TransientFilterType::Mininum;
  }

  TransientFilter filter(config);
  Progress progress;
  std::atomic_bool quit_progress(false);
  std::thread progress_thread([&progress, &quit_progress]() { runProrgess(progress, quit_progress); });

  filter.filter(cloud, &progress);
  const Cloud &transient = filter.transientResults();
  const Cloud &fixed = filter.fixedResults();

  quit_progress = true;
  progress_thread.join();
#else   // NEW_FILTER
  Cloud transient;
  Cloud fixed;
  cloud.findTransients(transient, fixed, merge_type, num_rays, colour);
#endif  // NEW_FILTER

  string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);

  transient.save(file_stub + "_transient.ply");
  fixed.save(file_stub + "_fixed.ply");
  return 0;
}
