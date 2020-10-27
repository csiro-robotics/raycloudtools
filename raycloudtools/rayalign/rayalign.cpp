// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raylib/raycloud.h"
#include "raylib/rayalignment.h"
#include "raylib/rayfinealignment.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayparse.h"

#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <complex>

void usage(int exit_code = 1)
{
  std::cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayalign raycloudA raycloudB" << std::endl;
  std::cout << "                             --nonrigid - nonrigid (quadratic) alignment" << std::endl;
  std::cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << std::endl;
  std::cout
    << "                             --local    - fine alignment only, assumes clouds are already approximately aligned"
    << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_a, cloud_b;
  ray::OptionalFlagArgument nonrigid("nonrigid", 'n'), is_verbose("verbose", 'v'), local("local", 'l');
  if (!ray::parseCommandLine(argc, argv, {&cloud_a, &cloud_b}, {&nonrigid, &is_verbose, &local}))
    usage();

  ray::Cloud clouds[2];
  if (!clouds[0].load(cloud_a.name()))
    usage();
  if (!clouds[1].load(cloud_b.name()))
    usage();

  bool local_only = local.isSet();
  bool non_rigid = nonrigid.isSet();
  bool verbose = is_verbose.isSet();
  if (verbose)
    ray::DebugDraw::init(argc, argv, "rayalign");
  if (!local_only)
  {
    alignCloud0ToCloud1(clouds, 0.5, verbose);
    if (verbose)
      clouds[0].save(cloud_a.nameStub() + "_coarse_aligned.ply");
  }

  ray::FineAlignment fineAlign(clouds, non_rigid, verbose);
  fineAlign.align();

  clouds[0].save(cloud_a.nameStub() + "_aligned.ply");
  return 0;
}
