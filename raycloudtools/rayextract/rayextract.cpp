// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayextractwoods.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(bool error=false)
{
  std::cout << "Extract feature into a text file structure" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayextract cloud.ply woods    - for low coverage scans (e.g. flyovers), where the trees are approximately inferred" << std::endl;
//  cout << " --extrapolate  - estimates tree distribution and adds trees where there is no evidence to the contrary" << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::DebugDraw::init(argc, argv, "rayextract");

  ray::KeyChoice extract_type({"woods"}); // "trees", "woods", "forest"});
  ray::FileArgument cloud_file; 
  ray::OptionalFlagArgument verbose("verbose", 'v');

  // three-way merge option
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &extract_type}, {&verbose}))
    usage();

  ray::Cloud cloud;
  cloud.load(cloud_file.name());

  if (extract_type.selectedKey() == "woods")
  {
    const double radius = 0.15;
    const double length = 1.0;
    ray::Wood woods(cloud, radius, length, verbose.isSet());
    woods.save(cloud_file.nameStub() + "_woods.txt");
  }
  return true;
}


