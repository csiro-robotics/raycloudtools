// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayparse.h"
#include "raylib/raycloudwriter.h"
#include "raylib/raycuboid.h"

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
  std::cout << "                  startpos 1,2,3         - splits based on start position, around plane 1,2,3" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays"
       << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  speed 1.0              - splits out rays when sensor moving above the given speed" << std::endl;
  std::cout << "                  time 1000 (or time 3 %)- splits at given time stamp (or percentage along)" << std::endl;
  std::cout << "                  box cx,cy,cz, rx,ry,rz - splits around an axis-aligned box at cx,cy,cz of given radii" << std::endl;
  std::cout << "                  grid wx,wy,wz          - splits into a 0,0,0 centred grid of files, cell width wx,wy,wz" << std::endl;
  exit(exit_code);
}

/// This is a helper function to aid in splitting the cloud while chunk-loading it. The purpose is to be able to
/// split clouds of any size, without running out of main memory. 
void split(const std::string &file_name, const std::string &in_name, 
           const std::string &out_name, std::function<bool(const ray::Cloud &cloud, int i)> is_outside)
{
  ray::Cloud cloud_buffer;
  ray::CloudWriter in_writer, out_writer;
  std::ofstream inside_ofs, outside_ofs;
  if (!in_writer.begin(in_name))
    usage();
  if (!out_writer.begin(out_name))
    usage();
  ray::Cloud in_chunk, out_chunk;

  /// move each ray into either the in_chunk or out_chunk, depending on the condition function is_outside
  auto per_chunk = [&cloud_buffer, &in_writer, &out_writer, &in_chunk, &out_chunk, &is_outside](
    std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, 
    std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in is_outside (by index). 
    cloud_buffer.starts = starts;
    cloud_buffer.ends = ends;
    cloud_buffer.times = times;
    cloud_buffer.colours = colours;

    for (int i = 0; i < (int)cloud_buffer.ends.size(); i++)
    {
      ray::Cloud &cloud = is_outside(cloud_buffer, i) ? out_chunk : in_chunk;
      cloud.addRay(cloud_buffer.starts[i], cloud_buffer.ends[i], cloud_buffer.times[i], cloud_buffer.colours[i]);
    }
    in_writer.writeChunk(in_chunk);
    out_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!ray::Cloud::read(file_name, per_chunk))
    usage(); 
  in_writer.end();
  out_writer.end();
}

/// Special case for splitting a box. 
void splitBox(const std::string &file_name, const std::string &in_name, const std::string &out_name, 
           const Eigen::Vector3d &centre, const Eigen::Vector3d &extents)
{
  ray::CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name))
    usage();
  if (!outside_writer.begin(out_name))
    usage();
  ray::Cloud in_chunk, out_chunk;

  /// This function assumes (correctly) that the passed in arguments are at end-of-life from readPly. 
  /// So it is valid to use std::move on them.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d p = ends[i] - centre;
      ray::Cuboid cuboid(-extents, extents);
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (cuboid.clipRay(start, end))
      {
        ray::RGBA col = colours[i];
        if (!cuboid.intersects(ends[i])) // mark as unbounded for the in_chunk
        {
          col.alpha = 0;
        }
        in_chunk.addRay(start, end, times[i], col);
        if (start != starts[i]) // start part is clipped
        {
          col.alpha = 0;
          out_chunk.addRay(starts[i], start, times[i], col);
        }
        if (ends[i] != end) // end part is clipped
        {
          out_chunk.addRay(end, ends[i], times[i], colours[i]);
        }
      }
      else
      {
        out_chunk.addRay(starts[i], ends[i], times[i], colours[i]);
      }
    }   
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!ray::readPly(file_name, true, per_chunk, 0))
    usage(); 

  inside_writer.end();
  outside_writer.end();
}

/// Special case for splitting based on a grid. 
void splitGrid(const std::string &file_name, const std::string &cloud_name, const Eigen::Vector3d &cell_width)
{
  ray::Cloud::CloudInfo info;
  ray::Cloud::getInfo(cloud_name, info);
  Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;
  
  Eigen::Vector3d minID(std::floor(min_bound[0]/cell_width[0]), std::floor(min_bound[1]/cell_width[1]), std::floor(min_bound[2]/cell_width[2]));
  Eigen::Vector3d maxID(std::ceil(max_bound[0]/cell_width[0]), std::ceil(max_bound[1]/cell_width[1]), std::ceil(max_bound[2]/cell_width[2]));
  Eigen::Vector3i minIndex = minID.cast<int>();
  Eigen::Vector3i maxIndex = maxID.cast<int>();
  Eigen::Vector3i dimensions = Eigen::Vector3i(1,1,1) + maxIndex - minIndex;
  int length = dimensions[0] * dimensions[1] * dimensions[2];
  if (length > 100000)
  {
    std::cout << "Error: grid is too many cells" << std::endl;
    usage();
  }

  std::vector<ray::CloudWriter> cells(length);
  std::vector<ray::Cloud> chunks(length);

  /// This function assumes (correctly) that the passed in arguments are at end-of-life from readPly. 
  /// So it is valid to use std::move on them.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    for (size_t i = 0; i < ends.size(); i++)
    {
      // how does the ray cross the different cells?
      // I guess I need to walk the cells (this again!!!)
      Eigen::Vector3d minI(std::floor(ends[i][0]/cell_width[0]), std::floor(ends[i][1]/cell_width[1]), std::floor(ends[i][2]/cell_width[2]));
      Eigen::Vector3d maxI(std::ceil(ends[i][0]/cell_width[0]), std::ceil(ends[i][1]/cell_width[1]), std::ceil(ends[i][2]/cell_width[2]));
      for (int x = minI[0]; x<=maxI[0]; x++)
      {
        for (int y = minI[1]; y<=maxI[1]; y++)
        {
          for (int z = minI[2]; z<=maxI[2]; z++)
          {
            int index = (x-minIndex[0]) + dimensions[0]*(y-minIndex[1]) + dimensions[0]*dimensions[1]*(z-minIndex[2]);
            if (cells[index].fileName() == "") // first time in this cell, so start writing to a new file
            {
              std::stringstream name;
              name << cloud_name << "_" << x << "_" << y << "_" << z << ".ply"; 
              cells[index].begin(name.str());
            } 
            // do actual clipping here.... 
            Eigen::Vector3d box_min((double)x*cell_width[0], (double)y*cell_width[1], (double)z*cell_width[2]);
            ray::Cuboid cuboid(box_min, box_min + cell_width);
            Eigen::Vector3d start = starts[i];
            Eigen::Vector3d end = ends[i];
            if (cuboid.clipRay(start, end))
            {
              ray::RGBA col = colours[i];
              if (!cuboid.intersects(ends[i])) // end point is outside, so mark an unbounded ray
              {
                col.alpha = 0;
              }
              chunks[index].addRay(start, end, times[i], col);
            }
          }
        }
      }
    }
    for (int i = 0; i<length; i++)
    {
      if (chunks[i].ends.size() > 0)
      {
        cells[i].writeChunk(chunks[i]);
        chunks[i].clear();
      }
    }       
  };
  if (!ray::Cloud::read(file_name, per_chunk))
    usage(); 

  for (int i = 0; i<length; i++)
  {
    if (cells[i].fileName() != "")
    {
      cells[i].end();
    }
  }  
}


// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument plane, colour(0.0, 1.0), startpos, raydir(-1.0, 1.0), box_centre, box_radius(0.0001, 1e10);
  ray::DoubleArgument time, alpha(0.0,1.0), range(0.0,1000.0), speed(0.0,1000.0);
  ray::KeyValueChoice choice({"plane", "time", "colour", "alpha", "startpos", "raydir", "range", "speed"}, 
                             {&plane,  &time,  &colour,  &alpha,  &startpos,  &raydir,  &range,  &speed});
  ray::FileArgument mesh_file;
  ray::TextArgument distance_text("distance"), time_text("time"), percent_text("%");
  ray::TextArgument box_text("box");
  ray::DoubleArgument mesh_offset;
  bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &choice});
  bool time_percent = ray::parseCommandLine(argc, argv, {&cloud_file, &time_text, &time, &percent_text});
  bool box_format = ray::parseCommandLine(argc, argv, {&cloud_file, &box_text, &box_centre, &box_radius});
  bool mesh_split = ray::parseCommandLine(argc, argv, {&cloud_file, &mesh_file, &distance_text, &mesh_offset});
  if (!standard_format && !box_format && !mesh_split && !time_percent)
  {
    usage();
  }

  const std::string in_name = cloud_file.nameStub() + "_inside.ply";
  const std::string out_name = cloud_file.nameStub() + "_outside.ply";
  const std::string rc_name = cloud_file.name(); // ray cloud name

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
    std::cout << "minimum time: " << min_time << " maximum time: " << max_time << ", difference: " 
              << max_time - min_time << std::endl;

    // now split based on this
    const double time_thresh = min_time + (max_time - min_time) * time.value()/100.0;
    split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
      return cloud.times[i] > time_thresh; 
    });
  }
  else if (box_format)
  {
    // Can't use cloud::split as sets are not mutually exclusive here.
    // we need to include rays that pass through the box. The intensity of these rays needs to be set to 0
    // so that they are treated as unbounded.
    splitBox(rc_name, in_name, out_name, box_centre.value(), box_radius.value());
  }  
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.times[i] > time.value(); 
      });
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.colours[i].alpha > c; 
      });
    }
    else if (parameter == "plane")
    {
      Eigen::Vector3d vec = plane.value() / plane.value().squaredNorm();
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.ends[i].dot(vec) > 1.0; 
      });
    }
    else if (parameter == "startpos")
    {
      Eigen::Vector3d vec = startpos.value() / startpos.value().squaredNorm();
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return cloud.starts[i].dot(vec) > 0.0; 
      });
    }
    else if (parameter == "raydir")
    {
      Eigen::Vector3d vec = raydir.value() / raydir.value().squaredNorm();
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 0.0;
      });
    }
    else if (parameter == "colour")
    {
      Eigen::Vector3d vec = colour.value() / colour.value().squaredNorm();
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                     (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 0.0;
      });
    }
    else if (parameter == "range")
    {
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool { 
        return (cloud.starts[i] - cloud.ends[i]).norm() > range.value(); 
      });
    }
    else if (parameter == "speed")
    {
      split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        if (i == 0)
        {
          return false;
        }
        return (cloud.starts[i] - cloud.starts[i - 1]).norm() / (cloud.times[i] - cloud.times[i - 1]) > speed.value();
      });
    }
  }
  return 0;
}
