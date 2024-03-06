// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "raylib/raycloud.h"
#include "raylib/raylaz.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raytrajectory.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Import a point cloud and trajectory file into a ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayimport pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << std::endl;
  std::cout << "                                           trajectoryfile is a text file using 'time x y z' format per line" << std::endl;
  std::cout << "rayimport pointcloudfile 0,0,0           - use 0,0,0 as the sensor location" << std::endl;
  std::cout << "rayimport pointcloudfile ray 0,0,-10     - use 0,0,-10 as the constant ray vector from start to point" << std::endl;
  std::cout << "                                        --max_intensity 100 - specify maximum intensity value (default 100)." << std::endl;
  std::cout << "                                                              0 sets all to full intensity (bounded rays)." << std::endl;
  std::cout << "                                        --remove_start_pos  - translate so first point is at 0,0,0" << std::endl;
  std::cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayImport(int argc, char *argv[])
{
  ray::DoubleArgument max_intensity(0.0, 10000, 100.0);
  ray::Vector3dArgument position, ray_vec;
  ray::TextArgument ray_text("ray");
  ray::OptionalKeyValueArgument max_intensity_option("max_intensity", 'm', &max_intensity);
  ray::OptionalFlagArgument remove("remove_start_pos", 'r');
  ray::FileArgument cloud_file, trajectory_file;
  bool standard_format =
    ray::parseCommandLine(argc, argv, { &cloud_file, &trajectory_file }, { &max_intensity_option, &remove });
  bool position_format =
    ray::parseCommandLine(argc, argv, { &cloud_file, &position }, { &max_intensity_option, &remove });
  bool ray_format =
    ray::parseCommandLine(argc, argv, { &cloud_file, &ray_text, &ray_vec }, { &max_intensity_option, &remove });
  if (!standard_format && !position_format && !ray_format)
    usage();

  if (ray_format && ray_vec.value().norm() == 0.0)
  {
    std::cerr << "Error: some ray cloud functions require rays to have a length. Please enter a non-zero vector for ray argument" << std::endl;
    usage();
  }
  ray::Cloud cloud;
  const std::string &traj_file = trajectory_file.name();
  // Sensors we use have 0 to 100 for normal output, and to 255 for special reflective surfaces
  double maximum_intensity = max_intensity.value();

  // load the trajectory first, it should fit into main memory
  ray::Trajectory trajectory;
  if (standard_format)
  {
    const std::string traj_end = traj_file.substr(traj_file.size() - 4);
    // allow the trajectory file to be in multiple different formats
    if (traj_end == ".ply" || traj_end == ".las" || traj_end == ".laz")
    {
      std::vector<Eigen::Vector3d> starts;
      std::vector<Eigen::Vector3d> ends;
      std::vector<double> times;
      std::vector<ray::RGBA> colours;
      if (traj_end == ".ply")
      {
        if (!ray::readPly(traj_file, starts, ends, times, colours, false))
          return false;
      }
      else
      {
        if (!ray::readLas(traj_file, ends, times, colours, maximum_intensity))
          return false;
      }
      trajectory.points() = std::move(ends);
      trajectory.times() = std::move(times);
    }
    else if (!trajectory.load(traj_file))
      usage();
  }

  std::string save_file = cloud_file.nameStub();
  if (cloud_file.nameExt() == "ply")
    save_file += "_raycloud";
  size_t num_bounded;
  std::ofstream ofs;
  ray::RayPlyBuffer buffer;
  if (!ray::writeRayCloudChunkStart(save_file + ".ply", ofs))
    usage();
  Eigen::Vector3d start_pos(0, 0, 0);
  bool has_warned = false;
  double min_time = std::numeric_limits<double>::max();
  double max_time = std::numeric_limits<double>::lowest();
  auto add_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<ray::RGBA> &colours) 
  {
    if (start_pos.squaredNorm() == 0.0)
    {
      start_pos = ends[0];
    }
    // user provides a single sensor location (e.g. for static scanners)
    if (position_format)
    {
      starts = ends;
      Eigen::Vector3d pos = position.value();
      for (auto &start : starts) 
      {
        start = pos;
      }
    }
    // user provides a constant ray vector
    // e.g. for an overhead aerial scan, if no trajectory is available
    else if (ray_format)
    {
      starts = ends;
      Eigen::Vector3d offset = -ray_vec.value();
      for (auto &start : starts) 
      {
        start += offset;
      }
    }
    // otherwise, a trajectory has been passed in
    else
    {
      // find the corresponding sensor locations for each point in the cloud
      trajectory.calculateStartPoints(times, starts);
      for (size_t i = 0; i < colours.size(); i++)
      {
        min_time = std::min(min_time, times[i]);
        max_time = std::max(max_time, times[i]);
        if (colours[i].alpha == 0 && ends[i][2] < starts[i][2])  // a nonreturn, we need to remove downward ones
        {
          Eigen::Vector3d dir = (ends[i] - starts[i]).normalized();
          const double minimal_distance_for_nonreturns = 0.1;
          ends[i] = starts[i] + dir * minimal_distance_for_nonreturns;
        }
      }
    }
    // option to remove the start position, for data that is in a global frame
    // this is particularly useful if we are storing the ray cloud positions using floats
    if (remove.isSet())
    {
      for (auto &end : ends) 
      {
        end -= start_pos;
      }
      for (auto &start : starts) 
      {
        start -= start_pos;
      }
    }
    if (maximum_intensity == 0.0)
    {
      for (auto &c : colours) 
      {
        c.alpha = 255;
      }
    }
    if (!ray::writeRayCloudChunk(ofs, buffer, starts, ends, times, colours, has_warned))
    {
      usage();
    }
  };
  Eigen::Vector3d *offset = remove.isSet() ? &start_pos : nullptr;
  if (cloud_file.nameExt() == "ply")
  {
    bool can_times_be_missing = position_format || ray_format;
    if (!ray::readPly(cloud_file.name(), false, add_chunk,
                      maximum_intensity, can_times_be_missing))  // special case of reading a non-ray-cloud ply
    {
      usage();
    }
  }
  else if (cloud_file.nameExt() == "laz" || cloud_file.nameExt() == "las")
  {
    if (!ray::readLas(cloud_file.name(), add_chunk, num_bounded, maximum_intensity, offset))
    {
      usage();
    }
  }
  else
  {
    std::cout << "Error converting unknown type: " << cloud_file.name() << std::endl;
    usage();
  }
  if (standard_format)
  {
    const float grace_period = 30.0;
    if (trajectory.times()[0] < min_time - grace_period)
    {
      std::cout << "trajectory begins " << min_time - trajectory.times()[0] << " s before first point cloud time" << std::endl;
    }
    if (trajectory.times().back() > max_time + grace_period)
    {
      std::cout << "trajectory ends " << trajectory.times().back() - max_time << " s after last point cloud time" << std::endl;
    }
    if (min_time < trajectory.times()[0]-grace_period || max_time > trajectory.times().back()+grace_period
     || min_time > trajectory.times().back() || max_time < trajectory.times()[0])
    {
      std::cerr.precision(10);
      std::cerr << "Error: trajectory times " << trajectory.times()[0] << "-" << trajectory.times().back() << 
        " do not span the point cloud times " << min_time << "-" << max_time << std::endl;
      usage();
    }
  }
  if (num_bounded == 0 && maximum_intensity > 0)
  {
    std::cout << "warning: all laz file intensities are 0." << std::endl;
    std::cout << "If your sensor lacks intensity information, set them to full using:" << std::endl;
    std::cout << "rayimport <point cloud> <trajectory file> --max_intensity 0" << std::endl;
  }
  ray::writeRayCloudChunkEnd(ofs);
  // if we remove the start position, then it is useful to print this value that is removed
  // so that the user hasn't lost information
  if (remove.isSet())
  {
    std::cout << "start position: " << start_pos.transpose() << " removed from all points" << std::endl;
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayImport, argc, argv);
}