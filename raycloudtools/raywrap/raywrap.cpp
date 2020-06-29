// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayconcavehull.h"
#include "raylib/rayconvexhull.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayply.h"
#include "raylib/rayutils.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>

// FIXME: Windows compatibility
#include <getopt.h>

void usage(int exitCode = 0)
{
  std::cout << "Extracts the ground surface as a mesh." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards"
       << std::endl;
  std::cout << "                               the 1.0 is the maximum curvature to bend to" << std::endl;
  std::cout << "--full                       - the full (slower) method accounts for overhangs." << std::endl;
  exit(exitCode);
}

int main(int argc, char *argv[])
{
  if (argc < 4 || argc > 5)
    usage();

  ray::DebugDraw::init(argc, argv, "ConcaveHull");
  std::string file = argv[1];
  std::string type = argv[2];
  double maximum_curvature = std::stod(argv[3]);
  bool overhangs = false;
  if (argc == 5)
  {
    if (std::string(argv[4]) == "--full" or std::string(argv[4]) == "-f")
      overhangs = true;
    else
      usage();
  }

  ray::Cloud cloud;
  cloud.load(file);
  cloud.removeUnboundedRays();
  if (file.substr(file.length() - 4) == ".ply")
    file = file.substr(0, file.length() - 4);

  if (overhangs)
  {
    ray::ConcaveHull concave_hull(cloud.ends);
    if (type == "inwards")
      concave_hull.growInwards(maximum_curvature);
    else if (type == "outwards")
      concave_hull.growOutwards(maximum_curvature);
    else if (type == "upwards")
      concave_hull.growUpwards(maximum_curvature);
    else if (type == "downwards")
      concave_hull.growTopDown(maximum_curvature);
    else
      usage();

    writePlyMesh(file + "_mesh.ply", concave_hull.mesh(), true);
  }
  else
  {
    ray::ConvexHull convexHull(cloud.ends);
    if (type == "inwards")
      convexHull.growInwards(maximum_curvature);
    else if (type == "outwards")
      convexHull.growOutwards(maximum_curvature);
    else if (type == "upwards")
      convexHull.growUpwards(maximum_curvature);
    else if (type == "downwards")
      convexHull.growTopDown(maximum_curvature);
    else
      usage();

    writePlyMesh(file + "_mesh.ply", convexHull.mesh, true);
  }


  std::cout << "Completed, output: " << file << "_mesh.ply" << std::endl;
  return 1;
}
