// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Reapply changes to a spatially decimated cloud back onto the full resolution cloud, on a per-voxel basis." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrestore decimated_cloud 10 cm full_res_cloud- decimated_cloud is a modified 10 cm decimation of full_res_cloud" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc != 5)
    usage();

  if (std::string(argv[3]) != "cm")
    usage();
  std::string decimated_file = argv[1];
  ray::Cloud decimated_cloud;
  decimated_cloud.load(decimated_file);

  std::string full_file = argv[4];
  ray::Cloud full_cloud;
  full_cloud.load(full_file);

  double voxel_width = 0.01 * std::stod(argv[2]);

  std::cout << "decimating" << std::endl;
  ray::Cloud full_decimated = full_cloud;
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
  full_decimated.decimate(voxel_width, &voxel_set); 

  std::cout << "decimated size: " << full_decimated.ends.size() << " altered size: " << decimated_cloud.ends.size() << std::endl;

  // how do we do this?
  // 1. detect alignment
  if (full_decimated.ends.size() == decimated_cloud.ends.size()) // no points added or removed, so look for a rigid transformation
  {
    // should we assume the point ordering is the same? I think we have to, otherwise we'd need to do a rayalign,
    // or at least use the covariance matrix to get the rigid transformation (bad if spherical!)

    // we could do a covariance matrix, if the eigenvalues are too similar, then do it based on distance weighted points, and run again?

    // simplest for now is to assume that they correspond, so pick 4 points (far apart) and get the transform of best fit from that
    // hang on, if they are spatially re-ordered, we can still correspond by time....
    // order both by time, to get in the same ordering.... they are both decimated anyway. 

    int is[3] = {0, (int)full_decimated.ends.size()/3, 2*(int)full_decimated.ends.size()/3};
    Eigen::Vector3d fullPs[3];
    Eigen::Vector3d decPs[3];
    Eigen::Vector3d midFull(0,0,0), midDec(0,0,0);
    for (int i = 0; i<3; i++)
    {
      fullPs[i] = full_decimated.ends[is[i]];
      decPs[i] = decimated_cloud.ends[is[i]];
      std::cout << "full time " << i << ": " << full_decimated.times[is[i]] - 1.513e9 << ", dec time: " << decimated_cloud.times[is[i]] - 1.513e9 << std::endl;
      midFull += fullPs[i]/3.0;
      midDec += decPs[i]/3.0;
    }
    for (int i = 1; i<3; i++)
    {
      fullPs[i] -= fullPs[0];
      decPs[i] -= decPs[0];
      std::cout << "fullPs[i] size: " << fullPs[i].norm() << ", decPs[i] size: " << decPs[i].norm() << std::endl;
    }

    // how to get rotation from two triangles? do it in two stages:
    Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(fullPs[1], decPs[1]);
    Eigen::Vector3d p1 = decPs[1].cross(decPs[2]);
    Eigen::Vector3d p2 = (fullPs[1].cross(fullPs[2]));
    Eigen::Quaterniond quat2 = Eigen::Quaterniond::FromTwoVectors(quat * p2, p1);
    Eigen::Quaterniond rotation = quat2 * quat;

    Eigen::Vector3d translation = midDec - rotation * midFull;
    full_cloud.transform(ray::Pose(translation, rotation), decimated_cloud.times[is[0]] - full_decimated.times[is[0]]);
  }
  else // points have been added or removed, so apply these additions and removals to the whole set of voxel points
  {
    std::set<Eigen::Vector3i, ray::Vector3iLess> test_set;
    // first build test_set:
    std::vector<Eigen::Vector3d> &points = decimated_cloud.ends;
    for (int64_t i = (int64_t)points.size()-1; i >= 0; i--)
    {
      Eigen::Vector3i place(int(std::floor(points[i][0] / voxel_width)), int(std::floor(points[i][1] / voxel_width)),
                            int(std::floor(points[i][2] / voxel_width)));
      test_set.insert(place);
    }

    // then for every point, we look up whether it exists in test_set, and add it if it does
    std::vector<Eigen::Vector3d> &ps = full_cloud.ends;
    std::vector<unsigned int> indices;
    for (unsigned int i = 0; i<ps.size(); i++)
    {
      Eigen::Vector3i place(int(std::floor(ps[i][0] / voxel_width)), int(std::floor(ps[i][1] / voxel_width)),
                            int(std::floor(ps[i][2] / voxel_width)));
      if (test_set.find(place) != test_set.end())
        indices.push_back(i);
    }
    std::cout << "indices found: " << indices.size() << std::endl;
    ray::Cloud new_cloud;
    new_cloud.starts.reserve(indices.size());
    new_cloud.ends.reserve(indices.size());
    new_cloud.times.reserve(indices.size());
    new_cloud.colours.reserve(indices.size());
    for (auto &i: indices)
      new_cloud.addRay(full_cloud, i);
    
    #define SUPPORT_ADDED_RAYS
    #if defined SUPPORT_ADDED_RAYS
    // what about additional points? can we make it work with a raycombine??
    std::vector<Eigen::Vector3d> &qs = decimated_cloud.ends;
    std::vector<unsigned int> indices2;
    for (unsigned int i = 0; i<qs.size(); i++)
    {
      Eigen::Vector3i place(int(std::floor(qs[i][0] / voxel_width)), int(std::floor(qs[i][1] / voxel_width)),
                            int(std::floor(qs[i][2] / voxel_width)));
      if (voxel_set.find(place) == voxel_set.end())
        indices2.push_back(i);
    }
    std::cout << "added " << indices2.size() << " extra points that are in the altered cloud" << std::endl;
    new_cloud.starts.reserve(indices.size() + indices2.size());
    new_cloud.ends.reserve(indices.size() + indices2.size());
    new_cloud.times.reserve(indices.size() + indices2.size());
    new_cloud.colours.reserve(indices.size() + indices2.size());
    for (auto &i: indices2)
      new_cloud.addRay(decimated_cloud, i);
    #endif

    full_cloud = new_cloud;
  }
  // finally do the same for colours, but only if the colours have been edited
  {

  }

  std::string file_stub = full_file;
  if (full_file.substr(full_file.length() - 4) == ".ply")
    file_stub = full_file.substr(0, full_file.length() - 4);
  full_cloud.save(file_stub + "_restored.ply");
  return true;
}
