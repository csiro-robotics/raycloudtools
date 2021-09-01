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
#include "raylib/rayparse.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>

// FIXME: Windows compatibility
#include <getopt.h>

void usage(int exit_code = 1)
{
  std::cout << "Extracts the ground surface as a mesh." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards"
       << std::endl;
  std::cout << "                               the 1.0 is the maximum curvature to bend to" << std::endl;
  std::cout << "--full                       - the full (slower) method accounts for overhangs." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  std::vector<Eigen::Vector2d> targets;
  std::vector<Eigen::Vector2d> vels, accs;
  double real_s = 3.1;
  double real_d = 2.5;
  // generate
  for (int i = 0; i<10; i++)
  {
    Eigen::Vector2d target(ray::random(-1.0, 1.0), ray::random(-1.0, 1.0));
    Eigen::Vector2d vel(ray::random(-1.0, 1.0), ray::random(-1.0, 1.0));
    Eigen::Vector2d acc = real_s*target + real_d*vel + Eigen::Vector2d(ray::random(-0.01, 0.01), ray::random(-0.01, 0.01));

    targets.push_back(target);
    vels.push_back(vel);
    accs.push_back(acc);
  }

  // recover
  Eigen::MatrixXd A(2*targets.size(), 2); 
  Eigen::VectorXd b(2*targets.size());
  for (size_t i = 0; i<targets.size(); i++)
  {
    for (int j = 0; j<2; j++)
    {
      A.row(2*i + j) = Eigen::Vector2d(targets[i][j], vels[i][j]);
      b[2*i + j] = accs[i][j];
    }
  }
  Eigen::Vector2d x = (A.transpose()*A).ldlt().solve(A.transpose()*b); 

  std::cout << "estimated s: " << x[0] << " vs real s: " << real_s << std::endl;
  std::cout << "estimated d: " << x[1] << " vs real d: " << real_d << std::endl;




  exit(1);
  ray::FileArgument cloud_file;
  ray::KeyChoice direction({"upwards", "downwards", "inwards", "outwards"});
  ray::DoubleArgument curvature;
  ray::OptionalFlagArgument full("full", 'f');
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &direction, &curvature}, {&full}))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();
  cloud.removeUnboundedRays();

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
    writePlyMesh(cloud_file.nameStub() + "_mesh.ply", convex_hull.mesh(), true);
  }

  std::cout << "Completed, output: " << cloud_file.nameStub() << "_mesh.ply" << std::endl;
  return 0;
}
