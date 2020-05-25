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
  std::cout << "Combines multiple ray clouds. Clouds are not moved but rays are omitted in the combined cloud according to "
          "the merge type specified."
       << std::endl;
  std::cout << "Outputs the combined cloud and the residual cloud of differences." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycombine min raycloud1 raycloud2 ... raycloudN 20 rays - combines into one cloud with minimal objects at "
          "differences"
       << std::endl;
  std::cout << "                                                           20 is the number of pass through rays to define "
          "a difference"
       << std::endl;
  std::cout
    << "           max    - maximal objects included. This is a form of volume intersection (rather than min: union)."
    << std::endl;
  std::cout << "           oldest - keeps the oldest geometry when there is a difference in later ray clouds." << std::endl;
  std::cout << "           newest - uses the newest geometry when there is a difference in newer ray clouds." << std::endl;
  std::cout << "           all    - combines as a simple concatenation, with all rays remaining (don't include 'xx rays')."
       << std::endl;
  std::cout << "raycombine basecloud min raycloud1 raycloud2 20 rays - 3-way merge, choses the changed geometry (from "
          "basecloud) at any differences. " << std::endl;
  std::cout << "For merge conflicts it uses the specified merge type." << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::DebugDraw::init(argc, argv, "raycombine");
  if (argc < 4)
    usage();
  bool threeway = std::string(argv[1])!="min" && std::string(argv[1])!="max" && std::string(argv[1])!="oldest" && 
    std::string(argv[1])!="newest" && std::string(argv[1])!="all";
  std::cout << "three way: " << threeway << std::endl;
  std::string merge_type = argv[threeway ? 2 : 1];
  bool concatenate = merge_type == "all";
  double num_rays = 0;
  if (concatenate)
  {
    if (std::string(argv[argc - 1]) == "rays")
      usage();
  }
  else
  {
    if (argc < 6 || std::string(argv[argc - 1]) != "rays")
      usage();
    num_rays = std::stod(argv[argc - 2]);
    argc -= 2;
  }

  int minI = threeway ? 3 : 2;
  if (threeway && argc != 5)
    usage();

  std::vector<std::string> files;
  for (int i = minI; i < argc; i++)
  {
    files.push_back(std::string(argv[i]));
    std::ifstream f(files.back().c_str());
    if (!f.good())
    {
      std::cout << "could not open file: " << files.back() << std::endl;
      usage();
    }
  }

  std::string file_stub = files[0];
  if (file_stub.substr(file_stub.length() - 4) == ".ply")
    file_stub = file_stub.substr(0, file_stub.length() - 4);

  std::vector<ray::Cloud> clouds(files.size());
  for (int i = 0; i < (int)files.size(); i++) clouds[i].load(files[i]);

  ray::Merger merger;
  if (threeway)
  {
    ray::Cloud base_cloud;
    base_cloud.load(argv[1]);
    merger.mergeThreeWay(base_cloud, clouds[0], clouds[1], merge_type, num_rays);
  }
  else if (concatenate)
  {
    ray::Cloud &result = merger.fixedCloud();
    for (auto &cloud : clouds)
    {
      result.starts.insert(result.starts.end(), cloud.starts.begin(), cloud.starts.end());
      result.ends.insert(result.ends.end(), cloud.ends.begin(), cloud.ends.end());
      result.times.insert(result.times.end(), cloud.times.begin(), cloud.times.end());
      result.colours.insert(result.colours.end(), cloud.colours.begin(), cloud.colours.end());
    }
  }
  else
  {
    merger.mergeMultipleClouds(clouds, merge_type, num_rays);
    merger.differenceCloud().save(file_stub + "_differences.ply");
  }
  merger.fixedCloud().save(file_stub + "_combined.ply");
  return true;
}
