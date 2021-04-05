// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayparse.h"
#include "raylib/raysplitter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits>

void usage(int exit_code = 1)
{
  std::cout << "Split a ray cloud relative to the supplied triangle mesh, generating two cropped ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysplit raycloud plane 10,0,0           - splits around plane at 10 m along x axis" << std::endl;
  std::cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << std::endl;
  std::cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << std::endl;
  std::cout << "                  meshfile distance 0.2  - splits raycloud at 0.2m from the meshfile surface" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays" << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  time 1000 (or time 3 %)- splits at given time stamp (or percentage along)" << std::endl;
  std::cout << "                  box rx,ry,rz           - splits around a centred axis-aligned box of the given radii" << std::endl;
  std::cout << "                  grid wx,wy,wz          - splits into a 0,0,0 centred grid of files, cell width wx,wy,wz" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  double max_val = std::numeric_limits<double>::max();
  ray::Vector3dArgument plane, colour(0.0, 1.0), raydir(-1.0, 1.0), box_radius(0.0001, max_val), cell_width(1.0, max_val);
  ray::DoubleArgument time, alpha(0.0,1.0), range(0.0,1000.0);
  ray::KeyValueChoice choice({"plane", "time", "colour", "alpha", "raydir", "range"}, 
                             {&plane,  &time,  &colour,  &alpha,  &raydir,  &range});
  ray::FileArgument mesh_file;
  ray::TextArgument distance_text("distance"), time_text("time"), percent_text("%");
  ray::TextArgument box_text("box"), grid_text("grid");
  ray::DoubleArgument mesh_offset;
  bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &choice});
  bool time_percent = ray::parseCommandLine(argc, argv, {&cloud_file, &time_text, &time, &percent_text});
  bool box_format = ray::parseCommandLine(argc, argv, {&cloud_file, &box_text, &box_radius});
  bool grid_format = ray::parseCommandLine(argc, argv, {&cloud_file, &grid_text, &cell_width});
  bool mesh_split = ray::parseCommandLine(argc, argv, {&cloud_file, &mesh_file, &distance_text, &mesh_offset});
  if (!standard_format && !box_format && !grid_format && !mesh_split && !time_percent)
  {
    usage();
  }

  const std::string in_name = cloud_file.nameStub() + "_inside.ply";
  const std::string out_name = cloud_file.nameStub() + "_outside.ply";
  const std::string rc_name = cloud_file.name(); // ray cloud name
  bool res = true;

  if (mesh_split) // I can't chunk load this one, so it will need to fit in RAM
  {
    ray::Cloud cloud; // used as a buffer when chunk loading
    if (!cloud.load(rc_name))
    {
      usage();
    }
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
    if (!ray::Cloud::read(cloud_file.name(), time_bounds))
      usage();
    std::cout << "Splitting cloud at " << (max_time - min_time) * time.value()/100.0 << 
      " seconds into the " << max_time - min_time << " time period of this ray cloud." << std::endl;

    // now split based on this
    const double time_thresh = min_time + (max_time - min_time) * time.value()/100.0;
    res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
      return cloud.times[i] > time_thresh; 
    });
  }
  else if (box_format)
  {
    // Can't use cloud::split as sets are not mutually exclusive here.
    // we need to include rays that pass through the box. The intensity of these rays needs to be set to 0
    // so that they are treated as unbounded.
    res = ray::splitBox(rc_name, in_name, out_name, Eigen::Vector3d(0,0,0), box_radius.value());
  }
  else if (grid_format)
  {
    // Can't use cloud::split as sets are not mutually exclusive here.
    // we need to include rays that pass through the box. The intensity of these rays needs to be set to 0
    // so that they are treated as unbounded.
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width.value());
  }    
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.times[i] > time.value(); 
      });
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.colours[i].alpha > c; 
      });
    }
    else if (parameter == "plane")
    {
      ray::splitPlane(rc_name, in_name, out_name, plane.value());
    }
    else if (parameter == "raydir")
    {
      Eigen::Vector3d vec = raydir.value() / raydir.value().squaredNorm();
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 1.0;
      });
    }
    else if (parameter == "colour")
    {
      Eigen::Vector3d vec = colour.value() / colour.value().squaredNorm();
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                     (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 1.0;
      });
    }
    else if (parameter == "range")
    {
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return (cloud.starts[i] - cloud.ends[i]).norm() > range.value(); 
      });
    }
  }
  if (!res)
    usage();
  return 0;
}
