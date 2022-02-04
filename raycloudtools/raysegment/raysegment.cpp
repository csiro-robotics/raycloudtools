// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/extraction/raysegment.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout << "Segment the ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysegment raycloud trees                  - colours each tree a different colour, assuming the trunk is visible" << std::endl;
  std::cout << "                        --max_diameter 0.9 - maximum trunk diameter for tree" << std::endl;
  std::cout << "                        --gradient 1.0     - maximum gradient for identifying ground points" << std::endl;
  std::cout << "                        --distance_limit 1 - maximum distance between points when building trees (beyond this are coloured black)" << std::endl;
  std::cout << "                        --height_min 2.0   - minimum height to count as a tree (under this are coloured black)" << std::endl;
  exit(exit_code);
}


// Colours the ray cloud based on the specified arguments
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::TextArgument trees_text("trees");
  ray::DoubleArgument max_diameter(0.01, 10.0), gradient(0.01, 20.0), distance_limit(0.01, 10.0), height_min(0.0, 100.0);
  ray::OptionalKeyValueArgument max_diameter_option("max_diameter", 'm', &max_diameter);
  ray::OptionalKeyValueArgument gradient_option("gradient", 'g', &gradient);
  ray::OptionalKeyValueArgument distance_limit_option("distance_limit", 'd', &distance_limit);
  ray::OptionalKeyValueArgument height_min_option("height_min", 'h', &height_min);
  const bool format = ray::parseCommandLine(argc, argv, {&cloud_file, &trees_text}, {&max_diameter_option, &gradient_option, &distance_limit_option, &height_min_option});
  if (!format)
    usage();
  
  std::string in_file = cloud_file.name();

  ray::Cloud cloud;
  if (!cloud.load(in_file))
    usage();
  double max_diam = max_diameter_option.isSet() ? max_diameter.value() : 0.9;
  double grad = gradient_option.isSet() ? gradient.value() : 1.0;
  double dist_limit = distance_limit_option.isSet() ? distance_limit.value() : 1.0;
  double height_minimum = height_min_option.isSet() ? height_min.value() : 2.0;
  segmentTrees(cloud, max_diam, grad, dist_limit, height_minimum);
  cloud.save(cloud_file.nameStub() + "_segmented.ply");
}
