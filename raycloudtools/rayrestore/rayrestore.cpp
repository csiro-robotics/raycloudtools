// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout << "Reapply changes to a decimated cloud back onto the full resolution cloud." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << " rayrestore decimated_cloud 10 cm full_cloud   - decimated_cloud is a modified 10 cm decimation of full_cloud" << std::endl;
  std::cout << " rayrestore decimated_cloud 10 rays full_cloud - decimated_cloud is a modified 'every tenth ray' decimation of full_cloud" << std::endl;
  std::cout << "Note: this tool does not work with raysmooth or temporal translations." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file, full_cloud_file;
  ray::DoubleArgument vox_width(0.1, 100.0);
  ray::IntArgument num_rays(1, 100);
  ray::ValueKeyChoice quantity({&vox_width, &num_rays}, {"cm", "rays"});
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &quantity, &full_cloud_file}))
    usage();


  ray::Cloud decimated_cloud;
  if (!decimated_cloud.load(cloud_file.name()))
    usage();

  ray::Cloud full_decimated; // we need a decimated version of the full cloud, to compare to
  bool spatial_decimation = quantity.selectedKey() == "cm";
#if 0
  ray::Cloud full_cloud;
  if (!full_cloud.load(full_cloud_file.name()))
    usage();

  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
  double voxel_width = 0;
  int ray_step = 0;
  std::cout << "decimating..." << std::endl;
  if (spatial_decimation)
  {
    voxel_width = 0.01 * vox_width.value();
    full_decimated = full_cloud;
    full_decimated.decimate(voxel_width, voxel_set); 
  }
  else
  {
    ray_step = num_rays.value();
    for (size_t i = 0; i < full_cloud.ends.size(); i += ray_step)
      full_decimated.addRay(full_cloud, i);
  }
#else
  std::vector<int64_t> subsample;
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;  
  full_decimated.reserve(decimated_cloud.ends.size()); // good guess at memory required

  auto decimate_cm = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    double width = 0.01 * vox_width.value();
    subsample.clear();
    voxelSubsample(ends, width, subsample, voxel_set);
    for (auto &id: subsample)
      full_decimated.addRay(starts[id], ends[id], times[id], colours[id]);
  };
  auto decimate_rays = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    size_t decimation = (size_t)num_rays.value();
    for (size_t i = 0, c = 0; i < ends.size(); i += decimation, c++)
      full_decimated.addRay(starts[i], ends[i], times[i], colours[i]);
  };

  bool success;
  if (spatial_decimation)
    success = ray::readPly(full_cloud_file.name(), true, decimate_cm, 0);
  else
    success = ray::readPly(full_cloud_file.name(), true, decimate_rays, 0);
  if (!success)
    usage();


  // temporary!
  ray::Cloud full_cloud;
  if (!full_cloud.load(full_cloud_file.name()))
    usage();

#endif
  std::cout << "full cloud decimated size: " << full_decimated.ends.size() << " modified cloud size: " << decimated_cloud.ends.size() << std::endl;

  struct Node
  {
    size_t index;
    double time;
  };
  std::vector<Node> decimated_nodes(decimated_cloud.times.size());
  std::vector<Node> full_decimated_nodes(full_decimated.times.size());
  for (size_t i = 0; i<decimated_cloud.times.size(); i++)
  {
    decimated_nodes[i].index = i;
    decimated_nodes[i].time = decimated_cloud.times[i];
  } 
  std::sort(decimated_nodes.begin(), decimated_nodes.end(), [](const Node &n1, const Node &n2){ return n1.time < n2.time; });
  for (size_t i = 0; i<full_decimated.times.size(); i++)
  {
    full_decimated_nodes[i].index = i;
    full_decimated_nodes[i].time = full_decimated.times[i];
  } 
  std::sort(full_decimated_nodes.begin(), full_decimated_nodes.end(), [](const Node &n1, const Node &n2){ return n1.time < n2.time; });
  

  double voxel_width = 0.01 * vox_width.value();
  int ray_step = num_rays.value();
  std::cout << "finding matching points" << std::endl;
  std::vector<Eigen::Vector2i> pairs; // corresponding ray indices between the two decimated clouds
  pairs.reserve(std::min(full_decimated.ends.size(), decimated_cloud.ends.size()));
  std::vector<size_t> added_ray_indices; // new rays added to the decimated_cloud
  added_ray_indices.reserve(decimated_cloud.ends.size());
  int j = 0;
  const double time_eps = 1e-6; // small enough to account for 200,000 rays per second, 
                                // but large enough to ignore file format/compression errors
  int num_removed_rays = 0;
  for (size_t i = 0; i<full_decimated_nodes.size(); i++)
  {
    while (decimated_nodes[j].time < full_decimated_nodes[i].time-time_eps && j<(int)decimated_nodes.size()-1)
    {
      added_ray_indices.push_back(decimated_nodes[j].index);
      j++;
    }
    if (std::abs(decimated_nodes[j].time - full_decimated_nodes[i].time) <= time_eps)
    {
      pairs.push_back(Eigen::Vector2i(full_decimated_nodes[i].index, decimated_nodes[j].index));
      j++;
    }
    else
    {
      if (spatial_decimation)
      {
        Eigen::Vector3d ps = full_decimated.ends[full_decimated_nodes[i].index];
        Eigen::Vector3i place(int(std::floor(ps[0] / voxel_width)), int(std::floor(ps[1] / voxel_width)),
                              int(std::floor(ps[2] / voxel_width)));
        voxel_set.erase(place);
      }
      num_removed_rays++;
    }
  }
  // finish adding the additional points
  while (j<(int)decimated_nodes.size()-1)
  {
    added_ray_indices.push_back(decimated_nodes[j].index);
    j++;
  }
  std::cout << "number of matched pairs: " << pairs.size() << ", number of removed rays: " << 
    num_removed_rays << ", number added: " << added_ray_indices.size() << std::endl;
 
  std::cout << "looking for a Euclidean transformation" << std::endl;
  // pick i and j indices from three points at different time points in the ray cloud
  // these represent the three corresponding pairs required to find the rigid transformation between
  // full_decimated cloud and decimated_cloud
  int is[3] = {pairs[0][0], pairs[pairs.size()/3][0], pairs[2*pairs.size()/3][0]};
  int js[3] = {pairs[0][1], pairs[pairs.size()/3][1], pairs[2*pairs.size()/3][1]};
  Eigen::Vector3d full_ps[3];
  Eigen::Vector3d dec_ps[3];
  Eigen::Vector3d mid_full(0,0,0), mid_dec(0,0,0);
  for (int i = 0; i<3; i++)
  {
    full_ps[i] = full_decimated.ends[is[i]];
    dec_ps[i] = decimated_cloud.ends[js[i]];
    mid_full += full_ps[i]/3.0;
    mid_dec += dec_ps[i]/3.0;
  }
  for (int i = 1; i<3; i++)
  {
    const double non_rigid_threshold = 0.01;
    full_ps[i] -= full_ps[0];
    dec_ps[i] -= dec_ps[0];
    if (abs(full_ps[i].norm() - dec_ps[i].norm()) > non_rigid_threshold)
      std::cout << "warning, matched points aren't a similar distance apart: " << full_ps[i].norm() <<
      ", " << dec_ps[i].norm() << " a non-rigid transform may have been applied. Results will be approximate." << std::endl;
  }

  // how to get rotation from two triangles? do it in two stages:
  Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(full_ps[1], dec_ps[1]);
  Eigen::Vector3d p1 = dec_ps[1].cross(dec_ps[2]);
  Eigen::Vector3d p2 = (full_ps[1].cross(full_ps[2]));
  Eigen::Quaterniond quat2 = Eigen::Quaterniond::FromTwoVectors(quat * p2, p1);
  Eigen::Quaterniond rotation = quat2 * quat;
  Eigen::Vector3d translation = mid_dec - rotation * mid_full;
  ray::Pose transform(translation, rotation);

  // now apply the estimated transform
  double rot_mag_sqr = ray::sqr(rotation.x()) + ray::sqr(rotation.y()) + ray::sqr(rotation.z());
  const double rotation_changed_threshold = 1e-8;
  const double translation_changed_threshold = 1e-8;
  if (rot_mag_sqr > ray::sqr(rotation_changed_threshold) || translation.squaredNorm() > ray::sqr(translation_changed_threshold))
  {
    std::cout << "transformation detected" << std::endl;
    std::cout << "translation: " << translation.transpose() << ", rotation quat: " << rotation.w() <<
      ", " << rotation.x() << ", " << rotation.y() << ", " << rotation.z() << std::endl;
  }
  else
  {
    std::cout << "no detected transformation of cloud" << std::endl;
    transform.position.setZero();
    transform.rotation = Eigen::Quaterniond::Identity();
  }

#if 1  
  std::ofstream ofs;
  if (!ray::writePlyChunkStart(full_cloud_file.nameStub() + "_restored.ply", ofs))
    usage();

  // By maintaining these buffers below, we avoid almost all memory fragmentation  
  ray::RayPlyBuffer buffer;
  ray::Cloud chunk;

  auto transfer_cm = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    chunk.clear();
    for (size_t i = 0; i<ends.size(); i++)
    {
      Eigen::Vector3i place(int(std::floor(ends[i][0] / voxel_width)), int(std::floor(ends[i][1] / voxel_width)),
                            int(std::floor(ends[i][2] / voxel_width)));
      if (voxel_set.find(place) != voxel_set.end())
        chunk.addRay(transform * starts[i], transform * ends[i], times[i], colours[i]);
    }
    ray::writePlyChunk(ofs, buffer, chunk.starts, chunk.ends, chunk.times, chunk.colours);
  };
  auto transfer_rays = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    int pair_index = 0;
    chunk.clear();
    for (int i = 0; i<static_cast<int>(ends.size()); i++)
    {
      int closest_index = (i+ray_step/2)/ray_step;
      while (pairs[pair_index][0] < closest_index && pair_index < (int)pairs.size()-1)
        pair_index++;
      if (pairs[pair_index][0] == closest_index)
        chunk.addRay(transform * starts[i], transform * ends[i], times[i], colours[i]);
    }
    ray::writePlyChunk(ofs, buffer, chunk.starts, chunk.ends, chunk.times, chunk.colours);
  };
  if (spatial_decimation)
    success = ray::readPly(full_cloud_file.name(), true, transfer_cm, 0);
  else
    success = ray::readPly(full_cloud_file.name(), true, transfer_rays, 0);
  if (!success)
    usage();

  std::cout << "added " << added_ray_indices.size() << " extra points that are in the modified cloud" << std::endl;
  chunk.clear(); 
  chunk.reserve(added_ray_indices.size());
  for (auto &ray_index: added_ray_indices)
    chunk.addRay(decimated_cloud, static_cast<int>(ray_index));
  ray::writePlyChunk(ofs, buffer, chunk.starts, chunk.ends, chunk.times, chunk.colours);

  ray::writePlyChunkEnd(ofs);
#else

  // then for every point in the full cloud, we have to look up whether it exists in the decimated clouds, 
  // and include it if it does
  std::cout << "applying full set of points to those found in decimated cloud" << std::endl;
  std::vector<Eigen::Vector3d> &ps = full_cloud.ends; // just a short hand
  std::vector<unsigned int> full_cloud_indices;
  full_cloud_indices.reserve(full_cloud.ends.size());
  if (spatial_decimation)
  {
    for (int i = 0; i<(int)ps.size(); i++)
    {
      Eigen::Vector3i place(int(std::floor(ps[i][0] / voxel_width)), int(std::floor(ps[i][1] / voxel_width)),
                            int(std::floor(ps[i][2] / voxel_width)));
      if (voxel_set.find(place) != voxel_set.end())
        full_cloud_indices.push_back(i);
    }
  }
  else
  {
    int pair_index = 0;
    for (int i = 0; i<(int)ps.size(); i++)
    {
      int closest_index = (i+ray_step/2)/ray_step;
      while (pairs[pair_index][0] < closest_index && pair_index < (int)pairs.size()-1)
        pair_index++;
      if (pairs[pair_index][0] == closest_index)
        full_cloud_indices.push_back(i);
    }
  }
  ray::Cloud new_cloud;
  new_cloud.reserve(full_cloud_indices.size());
  for (auto &i: full_cloud_indices)
    new_cloud.addRay(full_cloud, i);
  full_cloud = new_cloud;
  full_cloud.transform(transform, 0);
  

  // now, add back in the additional rays that were added to the modified cloud
  std::cout << "added " << added_ray_indices.size() << " extra points that are in the modified cloud" << std::endl;
  full_cloud.reserve(full_cloud.ends.size() + added_ray_indices.size());
  for (auto &ray_index: added_ray_indices)
    full_cloud.addRay(decimated_cloud, static_cast<int>(ray_index));

  full_cloud.save(full_cloud_file.nameStub() + "_restored.ply");
#endif
  return 0;
}
