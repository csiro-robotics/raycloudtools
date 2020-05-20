// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymerger.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/raydebugdraw.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Splits a raycloud into the transient rays and the fixed part" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raytransients min raycloud 20 rays - splits out positive transients (objects that have since moved)."
       << std::endl;
  std::cout << "                                     20 is number of pass through rays to classify as transient." << std::endl;
  std::cout << "              max    - finds negative transients, such as a hallway exposed when a door opens." << std::endl;
  std::cout << "              oldest - keeps the oldest geometry when there is a difference over time." << std::endl;
  std::cout << "              newest - uses the newest geometry when there is a difference over time." << std::endl;
  std::cout << " --colour     - also colours the clouds, to help tweak numRays. red: opacity, green: pass throughs, blue: "
          "planarity." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::DebugDraw::init(argc, argv, "raytransients");
  if (argc != 5 && argc != 6)
    usage();

  bool colour = false;
  if (argc == 6)
  {
    if (std::string(argv[5]) != "--colour" && std::string(argv[5]) != "-c")
      usage();
    colour = true;
  }
  double num_rays = std::stod(argv[3]);
  std::string merge_type = argv[1];
  if (merge_type != "min" && merge_type != "max" && merge_type != "oldest" && merge_type != "newest")
    usage();
  std::string file = argv[2];
  ray::Cloud cloud;
  cloud.load(file);

  ray::Merger merger;
  merger.mergeSingleCloud(cloud, merge_type, num_rays, colour);

  std::string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);

  merger.differences.save(file_stub + "_transient.ply");
  merger.result.save(file_stub + "_fixed.ply");
  return true;
}
