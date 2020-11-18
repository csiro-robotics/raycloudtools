// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits>

void usage(int exit_code = 1)
{
  std::cout << "Split a ray cloud relative to the supplied triangle mesh, generating two cropped ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysplit raycloud pos 10,0,0             - splits along x axis" << std::endl;
  std::cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << std::endl;
  std::cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << std::endl;
  std::cout << "                  meshfile distance 0.2  - splits raycloud at 0.2m from the meshfile surface" << std::endl;
  std::cout << "                  startpos 1,2,3         - splits based on start position, around plane 1,2,3" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays"
       << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  speed 1.0              - splits out rays when sensor moving above the given speed" << std::endl;
  std::cout << "                  time 1000 (or time 3 %)- splits at given time stamp (or percentage along)" << std::endl;
  exit(exit_code);
}

/// This is a helper function to aid in splitting the cloud while chunk-loading it. The purpose is to be able to
/// split clouds of any size, without running out of main memory. 
void split(ray::Cloud &cloud_buffer, const std::string &file_name, std::string &in_name, std::string &out_name, 
           std::function<bool(int i)> fptr)
{
  std::ofstream inside_ofs, outside_ofs;
  if (!ray::writePlyChunkStart(in_name, inside_ofs))
    usage();
  if (!ray::writePlyChunkStart(out_name, outside_ofs))
    usage();
  ray::RayPlyBuffer buffer;
  ray::Cloud in_chunk, out_chunk;

  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    cloud_buffer.starts = starts;
    cloud_buffer.ends = ends;
    cloud_buffer.times = times;
    cloud_buffer.colours = colours;

    for (int i = 0; i < (int)cloud_buffer.ends.size(); i++)
    {
      ray::Cloud &cloud = fptr(i) ? out_chunk : in_chunk;
      cloud.addRay(cloud_buffer.starts[i], cloud_buffer.ends[i], cloud_buffer.times[i], cloud_buffer.colours[i]);
    }
    ray::writePlyChunk(inside_ofs, buffer, in_chunk.starts, in_chunk.ends, in_chunk.times, in_chunk.colours);
    ray::writePlyChunk(outside_ofs, buffer, out_chunk.starts, out_chunk.ends, out_chunk.times, out_chunk.colours);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!ray::readPly(file_name, true, per_chunk, 0))
    usage(); 

  ray::writePlyChunkEnd(inside_ofs);
  ray::writePlyChunkEnd(outside_ofs);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument pos, colour(0.0, 1.0), startpos, raydir(-1.0, 1.0);
  ray::DoubleArgument time, alpha(0.0,1.0), range(0.0,1000.0), speed(0.0,1000.0);
  ray::KeyValueChoice choice({"pos", "time", "colour", "alpha", "startpos", "raydir", "range", "speed"}, 
                             {&pos,  &time,  &colour,  &alpha,  &startpos,  &raydir,  &range,  &speed});
  ray::FileArgument mesh_file;
  ray::TextArgument distance_text("distance"), time_text("time"), percent_text("%");
  ray::DoubleArgument mesh_offset;
  bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &choice});
  bool time_percent = ray::parseCommandLine(argc, argv, {&cloud_file, &time_text, &time, &percent_text});
  bool mesh_split = ray::parseCommandLine(argc, argv, {&cloud_file, &mesh_file, &distance_text, &mesh_offset});
  if (!standard_format && !mesh_split && !time_percent)
    usage();

  std::string in_name = cloud_file.nameStub() + "_inside.ply";
  std::string out_name = cloud_file.nameStub() + "_outside.ply";
  std::string rc_name = cloud_file.name(); // ray cloud name
  ray::Cloud cloud; // used as a buffer when chunk loading

  if (mesh_split) // I can't chunk load this one, so it will need to fit in RAM
  {
    if (!cloud.load(rc_name))
      usage();
    ray::Mesh mesh;
    ray::readPlyMesh(mesh_file.name(), mesh);
    ray::Cloud inside, outside;
    mesh.splitCloud(cloud, mesh_offset.value(), inside, outside);
    inside.save(in_name);
    outside.save(out_name);
  }
  else if (time_percent)
  {
    // chunk load the file just to get the time bounds
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::lowest();
    auto time_bounds = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &, std::vector<double> &times, std::vector<ray::RGBA> &)
    {
      for (auto &time: times)
      {
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
      }
    };
    if (!ray::readPly(cloud_file.name(), true, time_bounds, 0))
      usage();
    std::cout << "minimum time: " << min_time << " maximum time: " << max_time << ", difference: " 
              << max_time - min_time << std::endl;

    // now split based on this
    double time_thresh = min_time + (max_time - min_time) * time.value()/100.0;
    split(cloud, rc_name, in_name, out_name, [&](int i) -> bool { return cloud.times[i] > time_thresh; });
  }
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      split(cloud, rc_name, in_name, out_name, [&](int i) { return cloud.times[i] > time.value(); });
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      split(cloud, rc_name, in_name, out_name, [&](int i) { return cloud.colours[i].alpha > c; });
    }
    else if (parameter == "pos")
    {
      Eigen::Vector3d vec = pos.value() / pos.value().squaredNorm();
      split(cloud, rc_name, in_name, out_name, [&](int i) { return cloud.ends[i].dot(vec) > 1.0; });
    }
    else if (parameter == "startpos")
    {
      Eigen::Vector3d vec = startpos.value() / startpos.value().squaredNorm();
      split(cloud, rc_name, in_name, out_name, [&](int i) { return cloud.starts[i].dot(vec) > 0.0; });
    }
    else if (parameter == "raydir")
    {
      Eigen::Vector3d vec = raydir.value() / raydir.value().squaredNorm();
      split(cloud, rc_name, in_name, out_name, [&](int i) {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 0.0;
      });
    }
    else if (parameter == "colour")
    {
      Eigen::Vector3d vec = colour.value() / colour.value().squaredNorm();
      split(cloud, rc_name, in_name, out_name, [&](int i) {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                     (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 0.0;
      });
    }
    else if (parameter == "range")
    {
      split(cloud, rc_name, in_name, out_name, [&](int i) { 
        return (cloud.starts[i] - cloud.ends[i]).norm() > range.value(); 
      });
    }
    else if (parameter == "speed")
    {
      split(cloud, rc_name, in_name, out_name, [&](int i) {
        if (i == 0)
          return false;
        return (cloud.starts[i] - cloud.starts[i - 1]).norm() / (cloud.times[i] - cloud.times[i - 1]) > speed.value();
      });
    }
  }
  return 0;
}
