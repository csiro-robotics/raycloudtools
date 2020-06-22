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

  std::cout << "decimated size: " << full_decimated.ends.size() << " changed cloud size: " << decimated_cloud.ends.size() << std::endl;

  // This code assumes chronological ordering. If we move away from this, then 
  // we just need to sort the clouds chronologically first
  std::cout << "finding matching points" << std::endl;
  std::vector<Eigen::Vector2i> pairs;
  pairs.reserve(std::min(full_decimated.ends.size(), decimated_cloud.ends.size()));
  std::vector<int> missings;
  missings.reserve(decimated_cloud.ends.size());
  int j = 0;
  const double eps = 1e-6; // small enough to account for 200,000 rays per second, 
                           // but large enough to ignore file format/compression errors
  int count = 0;
  for (int i = 0; i<(int)full_decimated.ends.size() && j<(int)decimated_cloud.times.size(); i++)
  {
    while (decimated_cloud.times[j] < full_decimated.times[i]-eps && j<(int)decimated_cloud.times.size()-1)
    {
      missings.push_back(j);
      j++;
    }
    if (decimated_cloud.times[j] < full_decimated.times[i] + eps)
    {
      pairs.push_back(Eigen::Vector2i(i,j));
      j++;
    }
    else
    {
      Eigen::Vector3d ps = full_decimated.ends[i];
      Eigen::Vector3i place(int(std::floor(ps[0] / voxel_width)), int(std::floor(ps[1] / voxel_width)),
                            int(std::floor(ps[2] / voxel_width)));
      voxel_set.erase(place);
      count++;
    }
  }
  std::cout << "number of matched pairs: " << pairs.size() << ", number of removed points: " << count << ", number added: " << missings.size() << std::endl;
 
  std::cout << "looking for a Euclidean transformation" << std::endl;
  int is[3] = {pairs[0][0], pairs[pairs.size()/3][0], pairs[2*pairs.size()/3][0]};
  int js[3] = {pairs[0][1], pairs[pairs.size()/3][1], pairs[2*pairs.size()/3][1]};
  Eigen::Vector3d fullPs[3];
  Eigen::Vector3d decPs[3];
  Eigen::Vector3d midFull(0,0,0), midDec(0,0,0);
  for (int i = 0; i<3; i++)
  {
    fullPs[i] = full_decimated.ends[is[i]];
    decPs[i] = decimated_cloud.ends[js[i]];
    midFull += fullPs[i]/3.0;
    midDec += decPs[i]/3.0;
  }
  for (int i = 1; i<3; i++)
  {
    fullPs[i] -= fullPs[0];
    decPs[i] -= decPs[0];
    if ((fullPs[i] - decPs[i]).norm() > 0.01)
      std::cout << "warning, matched points aren't a very similar distance apart: " << fullPs[i].norm() <<
      ", " << decPs[i].norm() << " a non-rigid transform may have been applied" << std::endl;
  }

  // how to get rotation from two triangles? do it in two stages:
  Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(fullPs[1], decPs[1]);
  Eigen::Vector3d p1 = decPs[1].cross(decPs[2]);
  Eigen::Vector3d p2 = (fullPs[1].cross(fullPs[2]));
  Eigen::Quaterniond quat2 = Eigen::Quaterniond::FromTwoVectors(quat * p2, p1);
  Eigen::Quaterniond rotation = quat2 * quat;
  Eigen::Vector3d translation = midDec - rotation * midFull;

  // then for every point, we look up whether it exists in test_set, and add it if it does
  std::cout << "applying full set of points to those found in decimated cloud" << std::endl;
  std::vector<Eigen::Vector3d> &ps = full_cloud.ends;
  std::vector<unsigned int> indices;
  indices.reserve(full_cloud.ends.size());
  for (unsigned int i = 0; i<ps.size(); i++)
  {
    Eigen::Vector3i place(int(std::floor(ps[i][0] / voxel_width)), int(std::floor(ps[i][1] / voxel_width)),
                          int(std::floor(ps[i][2] / voxel_width)));
    if (voxel_set.find(place) != voxel_set.end())
      indices.push_back(i);
  }
  ray::Cloud new_cloud;
  new_cloud.starts.reserve(indices.size());
  new_cloud.ends.reserve(indices.size());
  new_cloud.times.reserve(indices.size());
  new_cloud.colours.reserve(indices.size());
  for (auto &i: indices)
    new_cloud.addRay(full_cloud, i);
  
  full_cloud = new_cloud;

  double r = ray::sqr(rotation.x()) + ray::sqr(rotation.y()) + ray::sqr(rotation.z());
  if (r > 1e-16 || translation.squaredNorm() > 1e-16)
  {
    std::cout << "transformation detected" << std::endl;
    std::cout << "translation: " << translation.transpose() << ", rotation: " << rotation.z() << std::endl;
    full_cloud.transform(ray::Pose(translation, rotation), decimated_cloud.times[is[0]] - full_decimated.times[is[0]]);
  }
  else
  {
    std::cout << "no detected transformation of cloud" << std::endl;
  }

  // now, add back in the missing ones
  std::cout << "added " << missings.size() << " extra points that are in the altered cloud" << std::endl;
  full_cloud.starts.reserve(full_cloud.starts.size() + missings.size());
  full_cloud.ends.reserve(full_cloud.ends.size() + missings.size());
  full_cloud.times.reserve(full_cloud.times.size() + missings.size());
  full_cloud.colours.reserve(full_cloud.colours.size() + missings.size());
  for (auto &missing: missings)
    full_cloud.addRay(decimated_cloud, missing);

  std::string file_stub = full_file;
  if (full_file.substr(full_file.length() - 4) == ".ply")
    file_stub = full_file.substr(0, full_file.length() - 4);
  full_cloud.save(file_stub + "_restored.ply");
  return true;
}
