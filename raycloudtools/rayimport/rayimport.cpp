// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "raylib/raycloud.h"
#include "raylib/rayply.h"
#include "raylib/raylaz.h"
#include "raylib/rayparse.h"
#include "raylib/raytrajectory.h"

void usage(int exit_code = 1)
{
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
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::DoubleArgument max_intensity(0.0, 10000);
  ray::Vector3dArgument position, ray_vec;
  ray::TextArgument ray_text("ray");
  ray::OptionalKeyValueArgument max_intensity_option("max_intensity", 'm', &max_intensity);
  ray::OptionalFlagArgument remove("remove_start_pos", 'r');
  ray::FileArgument cloud_file, trajectory_file;
  bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &trajectory_file}, {&max_intensity_option, &remove});
  bool position_format = ray::parseCommandLine(argc, argv, {&cloud_file, &position}, {&max_intensity_option, &remove});
  bool ray_format      = ray::parseCommandLine(argc, argv, {&cloud_file, &ray_text, &ray_vec}, {&max_intensity_option, &remove});
  if (!standard_format && !position_format && !ray_format)
    usage();

  ray::Cloud cloud;
  const std::string &traj_file = trajectory_file.name();
  // Sensors we use have 0 to 100 for normal output, and to 255 for special reflective surfaces
  double maximum_intensity = max_intensity_option.isSet() ? max_intensity.value() : 100.0;

  // load the trajectory first, it should fit into main memory
  ray::Trajectory trajectory;
  if (standard_format)
  {
    const std::string traj_end = traj_file.substr(traj_file.size() - 4);
    if (traj_end == ".ply")
    {
      std::vector<Eigen::Vector3d> starts;
      std::vector<Eigen::Vector3d> ends;
      std::vector<double> times;
      std::vector<ray::RGBA> colours;
      if (!ray::readPly(traj_file, starts, ends, times, colours, false))
        return false;
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
  Eigen::Vector3d start_pos(0,0,0);
  auto add_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    if (start_pos.squaredNorm() == 0.0)
    {
      start_pos = ends[0];
    }
    if (position_format)
    {
      starts = ends;
      Eigen::Vector3d pos = position.value();
      for (auto &start: starts)
        start = pos;      
    }
    else if (ray_format)
    {
      starts = ends;
      Eigen::Vector3d offset = -ray_vec.value();
      for (auto &start: starts)
        start += offset;         
    }
    else
    {
      trajectory.calculateStartPoints(times, starts);
      for (size_t i = 0; i<colours.size(); i++)
      {
        if (colours[i].alpha == 0 && ends[i][2] < starts[i][2]) // a nonreturn, we need to remove downward ones
        {
          Eigen::Vector3d dir = (ends[i] - starts[i]).normalized();
          const double minimal_distance_for_nonreturns = 0.1;
          ends[i] = starts[i] + dir*minimal_distance_for_nonreturns; 
        }
      }
    }
    if (remove.isSet())
    {
      for (auto &end: ends)
        end -= start_pos;
      for (auto &start: starts)
        start -= start_pos;
    }
    if (maximum_intensity == 0.0)
    {
      for (auto &c: colours)
        c.alpha = 255;
    }
    ray::writeRayCloudChunk(ofs, buffer, starts, ends, times, colours);
  };
  if (cloud_file.nameExt() == "ply")
  {
    if (!ray::readPly(cloud_file.name(), false, add_chunk, maximum_intensity)) // special case of reading a non-ray-cloud ply
      usage();
  }
  else if (cloud_file.nameExt() == "laz" || cloud_file.nameExt() == "las")
  {
    if (!ray::readLas(cloud_file.name(), add_chunk, num_bounded, maximum_intensity))
      usage();
  }
  else
  {
    std::cout << "Error converting unknown type: " << cloud_file.name() << std::endl;
    usage();
  }
  if (num_bounded == 0 && maximum_intensity > 0)
  {
    std::cout << "warning: all laz file intensities are 0." << std::endl;
    std::cout << "If your sensor lacks intensity information, set them to full using:" << std::endl;
    std::cout << "rayimport <point cloud> <trajectory file> --max_intensity 0" << std::endl;
  }  
  ray::writeRayCloudChunkEnd(ofs);
  if (remove.isSet())
  {
    std::cout << "start position: " << start_pos.transpose() << " removed from all points" << std::endl;
  }
  return 0;
}
