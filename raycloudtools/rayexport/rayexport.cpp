// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include "raylib/raycloud.h"
#include "raylib/raylaz.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raytrajectory.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Export a ray cloud into a point cloud amd trajectory file" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayexport raycloudfile.ply pointcloud.ply/.laz/.txt/.xyz trajectoryfile.ply/.txt - output in the chosen point cloud and trajectory formats" << std::endl;
  std::cout << "                           --traj_delta 0.1 - trajectory temporal decimation period in s. Default is 0.1" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayExport(int argc, char *argv[])
{
  ray::FileArgument raycloud_file, pointcloud_file, trajectory_file;
  ray::DoubleArgument traj_delta(0.0, 10000, 0.1);
  ray::OptionalKeyValueArgument delta_option("traj_delta", 't', &traj_delta);
  if (!ray::parseCommandLine(argc, argv, { &raycloud_file, &pointcloud_file, &trajectory_file }, { &delta_option }))
    usage();

  // Saving to a cloud file is fairly simple, we use chunk reading and writing:
  if (pointcloud_file.nameExt() == "laz")
  {
    ray::LasWriter las_writer(pointcloud_file.name());
    auto add_chunk = [&las_writer](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                   std::vector<double> &times,
                                   std::vector<ray::RGBA> &colours) { las_writer.writeChunk(ends, times, colours); };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
      usage();
  }
  else if (pointcloud_file.nameExt() == "ply")
  {
    ray::PointPlyBuffer buffer;
    std::ofstream ofs;
    if (!ray::writePointCloudChunkStart(pointcloud_file.name(), ofs))
      usage();
    bool has_warned = false;
    auto add_chunk = [&ofs, &buffer, &has_warned](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                     std::vector<double> &times, std::vector<ray::RGBA> &colours) {
      ray::writePointCloudChunk(ofs, buffer, ends, times, colours, has_warned);
    };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
      usage();
    ray::writePointCloudChunkEnd(ofs);
  }
  else if (pointcloud_file.nameExt() == "xyz" || pointcloud_file.nameExt() == "txt")
  {
    ray::PointPlyBuffer buffer;
    std::ofstream ofs;
    ofs.open(pointcloud_file.name(), std::ios::out);
    if (ofs.fail())
    {
      usage();
    }
    ofs << std::setprecision(4) << std::fixed;
    bool txt = pointcloud_file.nameExt() == "txt";
    if (txt)
    {
      ofs << "# point cloud text format. Comma delimited x,y,z,time,red,green,blue,alpha" << std::endl;
    }
    auto add_chunk = [&ofs, txt](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                            std::vector<double> &times, std::vector<ray::RGBA> &colours) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        if (txt)
        {
          ofs << ends[i][0] << "," << ends[i][1] << "," << ends[i][2] << "," << times[i] << "," << 
            (int)colours[i].red << "," << (int)colours[i].green << "," << (int)colours[i].blue << "," << (int)colours[i].alpha << std::endl;
        }
        else
        {
          ofs << ends[i][0] << " " << ends[i][1] << " " << ends[i][2] << " " << std::endl;
        }
      }    
    };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
    {
      usage();
    }
    ofs.close();
  }
  else
  {
    usage();
  }

  // saving the trajectory is more difficult. Firstly because we need to temporally decimate,
  // secondly because we need to sort the times, when saving to the txt file
  const double time_step = traj_delta.value();
  std::set<int64_t> time_slots;
  int64_t last_time_slot = std::numeric_limits<int64_t>::min();

  // if we are outputting to ply then we aren't sorting the times, just temporally decimating
  // that means we can still chunk-write the ply file, and the maximum memory is dictated by time_slots
  if (trajectory_file.nameExt() == "ply")
  {
    ray::PointPlyBuffer buffer;
    std::ofstream ofs;
    if (!ray::writePointCloudChunkStart(trajectory_file.name(), ofs))
      usage();
    ray::Cloud chunk;

    bool has_warned = false;
    auto decimate_time = [&time_slots, &ofs, &buffer, &chunk, &last_time_slot, time_step, &has_warned](
                           std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                           std::vector<double> &times, std::vector<ray::RGBA> &colours) {
      chunk.clear();
      for (size_t i = 0; i < ends.size(); i++)
      {
        const int64_t time_slot = static_cast<int64_t>(std::floor(times[i] / time_step));
        if (time_slot == last_time_slot)
          continue;
        if (time_slots.insert(time_slot).second)
        {
          chunk.starts.push_back(starts[i]);
          chunk.times.push_back(times[i]);
          chunk.colours.push_back(colours[i]);
        }
        last_time_slot = time_slot;
      }
      ray::writePointCloudChunk(ofs, buffer, chunk.starts, chunk.times, chunk.colours, has_warned);
    };
    if (!ray::readPly(raycloud_file.name(), true, decimate_time, 0))
      usage();
    ray::writePointCloudChunkEnd(ofs);
  }
  else if (trajectory_file.nameExt() == "txt")  // for text files we decimate and then sort
  {
    std::cout << "traj: " << trajectory_file.name() << std::endl;
    std::vector<ray::TrajectoryNode> traj_nodes;
    bool sorted = true;

    auto decimate_time = [&time_slots, &traj_nodes, &sorted, &last_time_slot, time_step](
                           std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                           std::vector<double> &times, std::vector<ray::RGBA> &) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        const int64_t time_slot = static_cast<int64_t>(std::floor(times[i] / time_step));
        if (time_slot == last_time_slot)
        {
          continue;
        }
        if (time_slots.insert(time_slot).second)
        {
          if (!traj_nodes.empty() && times[i] < traj_nodes.back().time)
          {
            sorted = false;
          }
          ray::TrajectoryNode traj_node;
          traj_node.time = times[i];
          traj_node.point = starts[i];
          traj_nodes.push_back(traj_node);
        }
        last_time_slot = time_slot;
      }
    };
    if (!ray::readPly(raycloud_file.name(), true, decimate_time, 0))
    {
      usage();
    }

    if (!sorted)
    {
      std::sort(traj_nodes.begin(), traj_nodes.end(),
                [](const ray::TrajectoryNode &a, const ray::TrajectoryNode &b) { return a.time < b.time; });
    }

    ray::saveTrajectory(traj_nodes, trajectory_file.name());
  }
  else
    usage();
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayExport, argc, argv);
}