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
  std::cout << "rayextract cloud.ply trunks    - estimate tree trunks and save to text file" << std::endl;
//  cout << " --extrapolate  - estimates tree distribution and adds trees where there is no evidence to the contrary" << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::DebugDraw::init(argc, argv, "rayextract");

  ray::KeyChoice extract_type({"trunks"}); // "tree", "trunks", "forest"});
  ray::FileArgument cloud_file; 
  ray::OptionalFlagArgument verbose("verbose", 'v');

  // three-way merge option
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &extract_type}, {&verbose}))
    usage();

  ray::Cloud cloud;
  cloud.load(cloud_file.name());

  if (extract_type.selectedKey() == "trunks")
  {
    const double radius = 0.15; // ~ /2 up to *2. So tree diameters 15 cm up to 60 cm 
    ray::Wood woods(cloud, radius, verbose.isSet());
    woods.save(cloud_file.nameStub() + "_trunks.txt");
  }
  return true;
}


