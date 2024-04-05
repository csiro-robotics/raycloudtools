// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayconcavehull.h"
#include "raylib/rayconvexhull.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/rayutils.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>

// FIXME: Windows compatibility

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Extracts the ground surface as a mesh." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards" << std::endl;
  std::cout << "                               the 1.0 is the maximum curvature to bend to" << std::endl;
  std::cout << "--full                       - the full (slower) method accounts for overhangs." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayWrap(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::KeyChoice direction({ "upwards", "downwards", "inwards", "outwards" });
  ray::DoubleArgument curvature;
  ray::OptionalFlagArgument full("full", 'f');
  if (!ray::parseCommandLine(argc, argv, { &cloud_file, &direction, &curvature }, { &full }))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();
  cloud.removeUnboundedRays();
  Eigen::Vector3d offset = cloud.removeStartPos();

  if (full.isSet())
  {
    ray::ConcaveHull concave_hull(cloud.ends);
    if (direction.selectedKey() == "inwards")
      concave_hull.growInwards(curvature.value());
    else if (direction.selectedKey() == "outwards")
      concave_hull.growOutwards(curvature.value());
    else if (direction.selectedKey() == "upwards")
      concave_hull.growUpwards(curvature.value());
    else if (direction.selectedKey() == "downwards")
      concave_hull.growDownwards(curvature.value());
    else
      usage();

    concave_hull.mesh().translate(offset);
    writePlyMesh(cloud_file.nameStub() + "_mesh.ply", concave_hull.mesh(), true);
  }
  else
  {
    ray::ConvexHull convex_hull(cloud.ends);
    if (direction.selectedKey() == "inwards")
      convex_hull.growInwards(curvature.value());
    else if (direction.selectedKey() == "outwards")
      convex_hull.growOutwards(curvature.value());
    else if (direction.selectedKey() == "upwards")
      convex_hull.growUpwards(curvature.value());
    else if (direction.selectedKey() == "downwards")
      convex_hull.growDownwards(curvature.value());
    else
      usage();

    convex_hull.mesh().reduce();
    convex_hull.mesh().translate(offset);
    writePlyMesh(cloud_file.nameStub() + "_mesh.ply", convex_hull.mesh(), true);
  }

  std::cout << "Completed, output: " << cloud_file.nameStub() << "_mesh.ply" << std::endl;
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayWrap, argc, argv);
}