// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raysplitter.h"
#include "raylib/rayforeststructure.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Split a ray cloud relative to the supplied geometry, generating two cropped ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysplit raycloud plane 10,0,0           - splits around plane at 10 m along x axis" << std::endl;
  std::cout << "                  colour                 - splits by colour, one cloud per colour" << std::endl;
  std::cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << std::endl;
  std::cout << "                  single_colour 255,0,0  - splits out a single colour, in 0-255 units" << std::endl;
  std::cout << "                  seg_colour             - splits to one cloud per colour, converting _segmented.ply colours to their index suffix" << std::endl;
  std::cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << std::endl;
  std::cout << "                  file distance 0.2      - splits raycloud at 0.2m from the (ply mesh or trees) file surface" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays" << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  time 1000 (or time 3 %)- splits at given time stamp (or percentage along)" << std::endl;
  std::cout << "                  box x,y,z rx,ry,rz     - splits around a given XYZ centred axis-aligned box of the given radii" << std::endl;  
  std::cout << "                  grid wx,wy,wz          - splits into a 0,0,0 centred grid of files, cell width wx,wy,wz. 0 for unused axes." << std::endl;
  std::cout << "                  grid wx,wy,wz 1        - same as above, but with a 1 metre overlap between cells." << std::endl;
  std::cout << "                  grid wx,wy,wz,wt       - splits into a grid of files, cell width wx,wy,wz and period wt. 0 for unused axes." << std::endl;
  std::cout << "                  capsule 1,2,3 10,11,12 5  - splits within a capsule using start, end and radius" << std::endl;
  // clang-format on
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int raySplit(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  double max_val = std::numeric_limits<double>::max();
  ray::Vector3dArgument plane, colour(0.0, 1.0), single_colour(0.0, 255.0), raydir(-1.0, 1.0),
    box_centre, box_radius(0.0, max_val), cell_width(0.0, max_val), capsule_start, capsule_end;
  ray::Vector4dArgument cell_width2(0.0, max_val);
  ray::DoubleArgument overlap(0.0, 10000.0);
  ray::DoubleArgument time, alpha(0.0, 1.0), range(0.0, 1000.0), capsule_radius(0.001, 1000.0);
  ray::KeyValueChoice choice({ "plane", "time", "colour", "single_colour", "alpha", "raydir", "range" },
                             { &plane, &time, &colour, &single_colour, &alpha, &raydir, &range });
  ray::FileArgument mesh_file, tree_file;
  ray::TextArgument distance_text("distance"), time_text("time"), percent_text("%");
  ray::TextArgument box_text("box"), grid_text("grid"), colour_text("colour"), seg_colour_text("seg_colour"), capsule_text("capsule");
  ray::DoubleArgument mesh_offset;
  bool standard_format = ray::parseCommandLine(argc, argv, { &cloud_file, &choice });
  bool colour_format = ray::parseCommandLine(argc, argv, { &cloud_file, &colour_text });
  bool seg_colour_format = ray::parseCommandLine(argc, argv, { &cloud_file, &seg_colour_text });
  bool time_percent = ray::parseCommandLine(argc, argv, { &cloud_file, &time_text, &time, &percent_text });
  bool box_format = ray::parseCommandLine(argc, argv, { &cloud_file, &box_text, &box_centre, &box_radius });
  bool grid_format = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width });
  bool grid_format2 = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width2 });
  bool grid_format3 = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width, &overlap });
  bool mesh_split = ray::parseCommandLine(argc, argv, { &cloud_file, &mesh_file, &distance_text, &mesh_offset });
  bool capsule_split =
    ray::parseCommandLine(argc, argv, { &cloud_file, &capsule_text, &capsule_start, &capsule_end, &capsule_radius });
  if (!standard_format && !colour_format && !seg_colour_format && !box_format && !grid_format && !grid_format2 && !grid_format3 &&
      !mesh_split && !time_percent && !capsule_split)
  {
    usage();
  }

  const std::string in_name = cloud_file.nameStub() + "_inside.ply";
  const std::string out_name = cloud_file.nameStub() + "_outside.ply";
  const std::string rc_name = cloud_file.name();  // ray cloud name
  bool res = true;

  // split the cloud around a capsule shape
  if (capsule_split)
  {
    res = ray::splitCapsule(rc_name, in_name, out_name, capsule_start.value(), capsule_end.value(), capsule_radius.value());
  }
  else if (colour_format)
  {
    res = ray::splitColour(cloud_file.name(), cloud_file.nameStub(), false);
  } 
  else if (seg_colour_format)
  {
    res = ray::splitColour(cloud_file.name(), cloud_file.nameStub(), true);
  }
  else if (mesh_split) 
  {
    if (mesh_file.nameExt() == "ply") // assume a mesh file
    {
      ray::Mesh mesh;
      ray::readPlyMesh(mesh_file.name(), mesh);
      if (!mesh.splitCloud(rc_name, mesh_offset.value(), in_name, out_name))
      {
        usage();
      }
    }
    else if (mesh_file.nameExt() == "txt") // assume a tree file
    {
      ray::ForestStructure forest;
      forest.load(mesh_file.name());
      ray::Cloud cloud;  // because forest splitCloud currently not chunk loaded
      if (!cloud.load(rc_name))
      {
        usage();
      }
      ray::Cloud inside, outside;
      forest.splitCloud(cloud, mesh_offset.value(), inside, outside);
      inside.save(in_name);
      outside.save(out_name);
    }
  }
  else if (time_percent)
  {
    // chunk load the file just to get the time bounds
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::lowest();
    auto time_bounds = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &, std::vector<double> &times,
                           std::vector<ray::RGBA> &) {
      for (auto &time : times)
      {
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
      }
    };
    if (!ray::Cloud::read(cloud_file.name(), time_bounds))
      usage();
    std::cout << "Splitting cloud at " << (max_time - min_time) * time.value() / 100.0 << " seconds into the "
              << max_time - min_time << " time period of this ray cloud." << std::endl;

    // now split based on this
    const double time_thresh = min_time + (max_time - min_time) * time.value() / 100.0;
    res = ray::split(rc_name, in_name, out_name,
                     [&](const ray::Cloud &cloud, int i) -> bool { return cloud.times[i] > time_thresh; });
  }
  else if (box_format)
  {
    Eigen::Vector3d extents = box_radius.value();
    for (int i = 0; i<3; i++) // use 0 for unbounded on an axis, for useability purposes, and to match the grid method
    {
      const double big_dimension = 1e7; // not too high just incase it causes precision issues inside clipRay
      if (extents[i] == 0.0)
      {
        extents[i] = big_dimension;
      }
    }
    res = ray::splitBox(rc_name, in_name, out_name, box_centre.value(), extents);
  }
  else if (grid_format)  // standard 3D grid of cuboids
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width.value());
  }
  else if (grid_format2)  // this is a 3+1D grid (space and time)
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width2.value());
  }
  else if (grid_format3)  // this is a 3D grid with a specified overlap
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width.value(), overlap.value());
  }
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      res = ray::split(rc_name, in_name, out_name,
                       [&](const ray::Cloud &cloud, int i) -> bool { return cloud.times[i] > time.value(); });
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      res = ray::split(rc_name, in_name, out_name,
                       [&](const ray::Cloud &cloud, int i) -> bool { return cloud.colours[i].alpha > c; });
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
    else if (parameter == "single_colour")  // split out a single colour
    {
      ray::RGBA col;
      col.red = (uint8_t)single_colour.value()[0];
      col.green = (uint8_t)single_colour.value()[1];
      col.blue = (uint8_t)single_colour.value()[2];
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        return !(cloud.colours[i].red == col.red && cloud.colours[i].green == col.green &&
                 cloud.colours[i].blue == col.blue);
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

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(raySplit, argc, argv);
}
