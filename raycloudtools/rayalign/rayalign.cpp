// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raylib/rayalignment.h"
#include "raylib/rayaxisalign.h"
#include "raylib/raycloud.h"
#include "raylib/rayfinealignment.h"
#include "raylib/rayparse.h"
#include "raylib/raypose.h"

#include <nabo/nabo.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << std::endl;
  std::cout << "This method is for when there is more than approximately 30% overlap between clouds." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayalign raycloudA raycloudB" << std::endl;
  std::cout << "                             --nonrigid - nonrigid (quadratic) alignment" << std::endl;
  std::cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << std::endl;
  std::cout << "                             --local    - fine alignment only, assumes clouds are already approximately aligned" << std::endl;
  std::cout << "rayalign raycloud  - axis aligns to the walls, placing the major walls at (0,0,0), biggest along y." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayAlign(int argc, char *argv[])
{
  ray::FileArgument cloud_a, cloud_b;
  ray::OptionalFlagArgument nonrigid("nonrigid", 'n'), is_verbose("verbose", 'v'), local("local", 'l');
  bool cross_align = ray::parseCommandLine(argc, argv, { &cloud_a, &cloud_b }, { &nonrigid, &is_verbose, &local });
  bool self_align = ray::parseCommandLine(argc, argv, { &cloud_a });
  if (!cross_align && !self_align)
    usage();

  std::string aligned_name = cloud_a.nameStub() + "_aligned.ply";
  if (self_align)
  {
    if (!ray::alignCloudToAxes(cloud_a.name(), aligned_name))
      usage();
  }
  else  // cross_align
  {
    ray::Pose transform;
    ray::Cloud clouds[2];
    if (!clouds[0].load(cloud_a.name()))
      usage();
    if (!clouds[1].load(cloud_b.name()))
      usage();

    // Here we pick two distant points in the cloud as an independent method of determining the total transformation
    // applied
    size_t min_i = 0, max_i = 0;
    for (size_t i = 0; i < clouds[0].ends.size(); i++)
    {
      if (clouds[0].ends[i][0] < clouds[0].ends[min_i][0])
        min_i = i;
      if (clouds[0].ends[i][0] > clouds[0].ends[max_i][0])
        max_i = i;
    }
    Eigen::Vector3d pos1 = clouds[0].ends[min_i];
    Eigen::Vector3d dir1 =
      Eigen::Vector3d(clouds[0].ends[max_i][0] - pos1[0], clouds[0].ends[max_i][1] - pos1[1], 0).normalized();

    bool local_only = local.isSet();
    bool non_rigid = nonrigid.isSet();
    bool verbose = is_verbose.isSet();
    if (!local_only)
    {
      alignCloud0ToCloud1(clouds, 0.5, verbose);
      if (verbose)
        clouds[0].save(cloud_a.nameStub() + "_coarse_aligned.ply");
    }

    ray::FineAlignment fineAlign(clouds, non_rigid, verbose);
    fineAlign.align();

    // Now we calculate the rigid transformation from the change in the position of the two points:
    Eigen::Vector3d pos2 = clouds[0].ends[min_i];
    Eigen::Vector3d dir2 =
      Eigen::Vector3d(clouds[0].ends[max_i][0] - pos2[0], clouds[0].ends[max_i][1] - pos2[1], 0).normalized();
    double angle = std::atan2((dir1.cross(dir2))[2], dir1.dot(dir2));
    Eigen::Vector3d rotated_pos1(pos1[0] * std::cos(angle) - pos1[1] * std::sin(angle),
                                 pos1[0] * std::sin(angle) + pos1[1] * std::cos(angle), pos1[2]);
    Eigen::Vector3d dif = pos2 - rotated_pos1;
    std::cout << "Transformation of " << cloud_a.nameStub() << ":" << std::endl;
    std::cout << "          rotation: (0, 0, " << angle * 180.0 / ray::kPi << ") degrees " << std::endl;
    std::cout << "  then translation: (" << dif.transpose() << ")" << std::endl;
    if (non_rigid)
    {
      std::cout << "This rigid transformation is approximate as a non-rigid transformation was applied" << std::endl;
    }

    clouds[0].save(aligned_name);
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayAlign, argc, argv);
}