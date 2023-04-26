// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/rayparse.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Reapply changes to a decimated cloud back onto the full resolution cloud." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << " rayrestore decimated_cloud 10 cm full_cloud   - decimated_cloud is a 10 cm decimation of full_cloud" << std::endl;
  std::cout << " rayrestore decimated_cloud 10 rays full_cloud - decimated_cloud is an 'every tenth ray' decimation of full_cloud" << std::endl;
  std::cout << "                         --restore_point_cloud - (-r) full_cloud is the original point cloud prior to rayimport. All ply fields retained" << std::endl;
  std::cout << "                         --point_cloud_start_removed - if using -r, was --remove_start_pos used in rayimport" << std::endl;
  std::cout << "Note: this tool does not work with raysmooth or temporal translations." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayRestore(int argc, char *argv[])
{
  ray::FileArgument cloud_file, full_cloud_file;
  ray::DoubleArgument vox_width(0.1, 100.0);
  ray::IntArgument num_rays(1, 100);
  ray::OptionalFlagArgument restore_point_cloud("restore_point_cloud", 'r');
  ray::OptionalFlagArgument start_removed("point_cloud_start_removed", 'p');
  ray::ValueKeyChoice quantity({ &vox_width, &num_rays }, { "cm", "rays" });
  if (!ray::parseCommandLine(argc, argv, { &cloud_file, &quantity, &full_cloud_file }, {&restore_point_cloud, &start_removed}))
    usage();
  const bool spatial_decimation = quantity.selectedKey() == "cm";
  const double voxel_width = 0.01 * vox_width.value();
  const int ray_step = num_rays.value();

  // This function uses chunk loading to avoid the full resolution cloud being in memory

  // Firstly, load the decimated cloud. We assume that this can fit in RAM
  ray::Cloud decimated_cloud;
  if (!decimated_cloud.load(cloud_file.name()))
    usage();

  ray::Cloud full_decimated;       // we need a decimated version of the full cloud, to compare to
  std::vector<int64_t> subsample;  // single buffer minimises memory allocations
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
  full_decimated.reserve(decimated_cloud.ends.size());  // good guess at memory required

  // decimation functions
  bool start_found = false;
  Eigen::Vector3d start_pos(0,0,0);
  auto decimate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    if (start_removed.isSet())
    {
      if (!start_found)
      {
        start_pos = ends[0];
        start_found = true;
        std::cout << "start position that was removed " << start_pos.transpose() << " has been restored." << std::endl; 
      }
      for (size_t i = 0; i<ends.size(); i++)
      {
        ends[i] -= start_pos; // remove start pos again so that we can pair up with the decimated cloud
        starts[i] -= start_pos;
      }
    }
    if (spatial_decimation)
    {
      subsample.clear();
      voxelSubsample(ends, voxel_width, subsample, voxel_set);
      
      for (auto &id : subsample) 
      {
        full_decimated.addRay(starts[id], ends[id], times[id], colours[id]);
      }
    }
    else
    {
      const int num_points = static_cast<int>(ends.size());
      for (int i = 0, c = 0; i < num_points; i += ray_step, c++)
        full_decimated.addRay(starts[i], ends[i], times[i], colours[i]);
    }
  };
  if (!ray::readPly(full_cloud_file.name(), !restore_point_cloud.isSet(), decimate, 0))
    usage();

  std::cout << "full cloud decimated size: " << full_decimated.ends.size()
            << " modified cloud size: " << decimated_cloud.ends.size() << std::endl;

  // Next we need to order the rays by time. Rather than shifting large data, we sort the indices:
  struct Node
  {
    size_t index;
    double time;
  };
  std::vector<Node> decimated_nodes(decimated_cloud.times.size());
  std::vector<Node> full_decimated_nodes(full_decimated.times.size());
  for (size_t i = 0; i < decimated_cloud.times.size(); i++)
  {
    decimated_nodes[i].index = i;
    decimated_nodes[i].time = decimated_cloud.times[i];
  }
  std::sort(decimated_nodes.begin(), decimated_nodes.end(),
            [](const Node &n1, const Node &n2) { return n1.time < n2.time; });
  for (size_t i = 0; i < full_decimated.times.size(); i++)
  {
    full_decimated_nodes[i].index = i;
    full_decimated_nodes[i].time = full_decimated.times[i];
  }
  std::sort(full_decimated_nodes.begin(), full_decimated_nodes.end(),
            [](const Node &n1, const Node &n2) { return n1.time < n2.time; });

  // Now find matching points by time. We assume that accurate time is a unique identifier per point
  std::cout << "finding matching points" << std::endl;
  std::vector<Eigen::Vector2i> pairs;  // corresponding ray indices between the two decimated clouds
  pairs.reserve(std::min(full_decimated.ends.size(), decimated_cloud.ends.size()));
  std::vector<size_t> added_ray_indices;  // new rays added to the decimated_cloud
  added_ray_indices.reserve(decimated_cloud.ends.size());
  int j = 0;
  const double time_eps = 1e-7;  // small enough to account for 200,000 rays per second,
                                 // but large enough to ignore file format/compression errors
  int num_removed_rays = 0;
  for (size_t i = 0; i < full_decimated_nodes.size(); i++)
  {
    // rays added into the decimated_cloud
    while (decimated_nodes[j].time < full_decimated_nodes[i].time - time_eps && j < (int)decimated_nodes.size() - 1)
    {
      added_ray_indices.push_back(decimated_nodes[j].index);
      j++;
    }
    // matching points
    if (std::abs(decimated_nodes[j].time - full_decimated_nodes[i].time) <= time_eps)
    {
      pairs.push_back(Eigen::Vector2i(full_decimated_nodes[i].index, decimated_nodes[j].index));
      j++;
    }
    // rays removed from the full_decimated cloud
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
  while (j < (int)decimated_nodes.size() - 1)
  {
    added_ray_indices.push_back(decimated_nodes[j].index);
    j++;
  }
  std::cout << "number of matched pairs: " << pairs.size() << ", number of removed rays: " << num_removed_rays
            << ", number added: " << added_ray_indices.size() << std::endl;

  // Now find the Euclidan transformation that transforms the point pairs
  std::cout << "looking for a Euclidean transformation" << std::endl;
  // pick i and j indices from three points at different time points in the ray cloud
  // these represent the three corresponding pairs required to find the rigid transformation between
  // full_decimated cloud and decimated_cloud
  ray::Pose transform;
  transform.position.setZero();
  transform.rotation = Eigen::Quaterniond::Identity();
  // only estimate a transform if there are a sufficient number of pairs
  if (pairs.size() >= 6)
  {
    const int is[3] = { pairs[0][0], pairs[pairs.size() / 3][0], pairs[2 * pairs.size() / 3][0] };
    const int js[3] = { pairs[0][1], pairs[pairs.size() / 3][1], pairs[2 * pairs.size() / 3][1] };
    Eigen::Vector3d full_ps[3];  // a triangle in the full_decimated cloud
    Eigen::Vector3d dec_ps[3];   // a triangle in the decimated_cloud
    Eigen::Vector3d mid_full(0, 0, 0), mid_dec(0, 0, 0);
    for (int i = 0; i < 3; i++)
    {
      full_ps[i] = full_decimated.ends[is[i]];
      dec_ps[i] = decimated_cloud.ends[js[i]];
      mid_full += full_ps[i] / 3.0;
      mid_dec += dec_ps[i] / 3.0;
    }
    for (int i = 1; i < 3; i++)
    {
      const double non_rigid_threshold = 0.01;
      full_ps[i] -= full_ps[0];
      dec_ps[i] -= dec_ps[0];
      if (std::abs(full_ps[i].norm() - dec_ps[i].norm()) > non_rigid_threshold)
      {
        std::cout << "warning, matched points aren't a similar distance apart: " << full_ps[i].norm() << ", "
                  << dec_ps[i].norm() << " a non-rigid transform may have been applied. Results will be approximate."
                  << std::endl;
      }
    }

    // how to get rotation from two triangles? do it in two stages:
    const Eigen::Quaterniond quat = Eigen::Quaterniond::FromTwoVectors(full_ps[1], dec_ps[1]);
    const Eigen::Vector3d normal1 = dec_ps[1].cross(dec_ps[2]);
    const Eigen::Vector3d normal2 = full_ps[1].cross(full_ps[2]);
    const Eigen::Quaterniond quat2 = Eigen::Quaterniond::FromTwoVectors(quat * normal2, normal1);
    const Eigen::Quaterniond rotation = quat2 * quat;
    const Eigen::Vector3d translation = mid_dec - rotation * mid_full;
    transform.position = translation;
    transform.rotation = rotation;

    // set transformation to identity, if it is very close. This makes the typical case more accurate
    const double rot_mag_sqr = ray::sqr(rotation.x()) + ray::sqr(rotation.y()) + ray::sqr(rotation.z());
    const double rotation_changed_threshold = 1e-8;
    const double translation_changed_threshold = 1e-8;
    if (rot_mag_sqr > ray::sqr(rotation_changed_threshold) ||
        translation.squaredNorm() > ray::sqr(translation_changed_threshold))
    {
      std::cout << "transformation detected" << std::endl;
      std::cout << "translation: " << translation.transpose() << ", rotation quat: " << rotation.w() << ", "
                << rotation.x() << ", " << rotation.y() << ", " << rotation.z() << std::endl;
    }
    else
    {
      std::cout << "no detected transformation of cloud" << std::endl;
    }
  }

  int pair_index = 0;
  size_t extra_rays = 0;
  auto convert = [&](int index, Eigen::Vector3d &start, Eigen::Vector3d &end, double &time, ray::RGBA &colour) -> int
  {
    if (extra_rays < added_ray_indices.size())
    {
      size_t id = added_ray_indices[extra_rays++];
      start = decimated_cloud.starts[id] + start_pos;
      end = decimated_cloud.ends[id] + start_pos;
      time = decimated_cloud.times[id];
      colour = decimated_cloud.colours[id];
      return 2;
    }
    if (spatial_decimation)
    {
      start -= start_pos;
      end -= start_pos;
      Eigen::Vector3i place(int(std::floor(end[0] / voxel_width)), int(std::floor(end[1] / voxel_width)),
                            int(std::floor(end[2] / voxel_width)));
      if (voxel_set.find(place) != voxel_set.end())
      {
        start = start_pos + transform * start;
        end = start_pos + transform * end;
        return 1;
      }
    }
    else
    {
      int closest_index = (index + ray_step / 2) / ray_step;
      while (pairs[pair_index][0] < closest_index && pair_index < (int)pairs.size() - 1) 
      {
        pair_index++;
      }
      if (pairs[pair_index][0] == closest_index)
      {
        start = start_pos + transform * (start - start_pos);
        end = start_pos + transform * (end - start_pos);
        return 1;
      }
    }
    return 0;
  };
  if (!ray::readGeneralPly(full_cloud_file.name(), full_cloud_file.nameStub() + "_restored.ply", !restore_point_cloud.isSet(), nullptr, convert, 0))
  {
    usage();
  }
  std::cout << "added " << added_ray_indices.size() << " extra points that are in the modified cloud" << std::endl;
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayRestore, argc, argv);
}
