// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <iostream>
#include <limits>
#include "raysplitter.h"
#include "raycloudwriter.h"
#include "raycuboid.h"

namespace ray
{
/// This is a helper function to aid in splitting the cloud while chunk-loading it. The purpose is to be able to
/// split clouds of any size, without running out of main memory. 
bool split(const std::string &file_name, const std::string &in_name, 
           const std::string &out_name, std::function<bool(const Cloud &cloud, int i)> is_outside)
{
  Cloud cloud_buffer;
  CloudWriter in_writer, out_writer;
  std::ofstream inside_ofs, outside_ofs;
  if (!in_writer.begin(in_name))
    return false;
  if (!out_writer.begin(out_name))
    return false;
  Cloud in_chunk, out_chunk;

  /// move each ray into either the in_chunk or out_chunk, depending on the condition function is_outside
  auto per_chunk = [&cloud_buffer, &in_writer, &out_writer, &in_chunk, &out_chunk, &is_outside](
    std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, 
    std::vector<double> &times, std::vector<RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in is_outside (by index). 
    cloud_buffer.starts = starts;
    cloud_buffer.ends = ends;
    cloud_buffer.times = times;
    cloud_buffer.colours = colours;

    for (int i = 0; i < (int)cloud_buffer.ends.size(); i++)
    {
      Cloud &cloud = is_outside(cloud_buffer, i) ? out_chunk : in_chunk;
      cloud.addRay(cloud_buffer.starts[i], cloud_buffer.ends[i], cloud_buffer.times[i], cloud_buffer.colours[i]);
    }
    in_writer.writeChunk(in_chunk);
    out_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!Cloud::read(file_name, per_chunk))
    return false; 
  in_writer.end();
  out_writer.end();
  return true;
}

/// Special case for splitting a plane. 
bool splitPlane(const std::string &file_name, const std::string &in_name, const std::string &out_name, const Eigen::Vector3d &plane)
{
  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name))
    return false;
  if (!outside_writer.begin(out_name))
    return false;
  Cloud in_chunk, out_chunk;

  /// This function assumes (correctly) that the passed in arguments are at end-of-life from readPly. 
  /// So it is valid to use std::move on them.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    Eigen::Vector3d p = plane / plane.dot(plane);
    for (size_t i = 0; i < ends.size(); i++)
    {
      double d1 = starts[i].dot(p) - 1.0;
      double d2 = ends[i].dot(p) - 1.0;
      if (d1*d2 > 0.0) // start and end are on the same side of the plane, so don't split...
      {
        Cloud &chunk = d1 > 0.0 ? out_chunk : in_chunk;
        chunk.addRay(starts[i], ends[i], times[i], colours[i]);
      }
      else // split the ray...
      {
        RGBA col = colours[i];
        col.alpha = 0;
        Eigen::Vector3d mid = starts[i] + (ends[i] - starts[i]) * d1/(d1-d2);
        if (d1 > 0.0)
        {
          out_chunk.addRay(starts[i], mid, times[i], col);
          in_chunk.addRay(mid, ends[i], times[i], colours[i]);
        }
        else
        {
          in_chunk.addRay(starts[i], mid, times[i], col);
          out_chunk.addRay(mid, ends[i], times[i], colours[i]);
        }
      }
    }   
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!readPly(file_name, true, per_chunk, 0))
    return false; 

  inside_writer.end();
  outside_writer.end();
  return true;
}

/// Special case for splitting a box. 
bool splitBox(const std::string &file_name, const std::string &in_name, const std::string &out_name, 
           const Eigen::Vector3d &centre, const Eigen::Vector3d &extents)
{
  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name))
    return false;
  if (!outside_writer.begin(out_name))
    return false;
  Cloud in_chunk, out_chunk;

  /// This function assumes (correctly) that the passed in arguments are at end-of-life from readPly. 
  /// So it is valid to use std::move on them.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    Cuboid cuboid(centre - extents, centre + extents);
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (cuboid.clipRay(start, end))
      {
        RGBA col = colours[i];
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
  if (!readPly(file_name, true, per_chunk, 0))
    return false; 

  inside_writer.end();
  outside_writer.end();
  return true;
}

/// Special case for splitting based on a grid. 
bool splitGrid(const std::string &file_name, const std::string &cloud_name_stub, const Eigen::Vector3d &cell_width)
{
  Cloud::CloudInfo info;
  Cloud::getInfo(cloud_name_stub + ".ply", info);
  Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;
  
  Eigen::Vector3d minID(std::floor(0.5 + min_bound[0]/cell_width[0]), std::floor(0.5 + min_bound[1]/cell_width[1]), std::floor(0.5 + min_bound[2]/cell_width[2]));
  Eigen::Vector3d maxID(std::ceil(0.5 + max_bound[0]/cell_width[0]), std::ceil(0.5 + max_bound[1]/cell_width[1]), std::ceil(0.5 + max_bound[2]/cell_width[2]));
  Eigen::Vector3i minIndex = minID.cast<int>();
  Eigen::Vector3i maxIndex = maxID.cast<int>();
  Eigen::Vector3i dimensions = maxIndex - minIndex;
  int length = dimensions[0] * dimensions[1] * dimensions[2];
  if (length > 1000000)
  {
    std::cout << "Error: grid is too many cells" << std::endl;
    return false;
  }

  std::vector<CloudWriter> cells(length);
  std::vector<Cloud> chunks(length);

  /// This function assumes (correctly) that the passed in arguments are at end-of-life from readPly. 
  /// So it is valid to use std::move on them.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours)
  {
    // I move these into the cloud buffer, so that they can be indexed easily in fptr (by index). 
    for (size_t i = 0; i < ends.size(); i++)
    {
      // how does the ray cross the different cells?
      // I guess I need to walk the cells (this again!!!)
      Eigen::Vector3d from = Eigen::Vector3d(0.5, 0.5, 0.5) + starts[i].cwiseQuotient(cell_width);
      Eigen::Vector3d to = Eigen::Vector3d(0.5, 0.5, 0.5) + ends[i].cwiseQuotient(cell_width);
      Eigen::Vector3d pos0 = minVector(from, to);
      Eigen::Vector3d pos1 = maxVector(from, to);
      Eigen::Vector3i minI = Eigen::Vector3d(std::floor(pos0[0]), std::floor(pos0[1]), std::floor(pos0[2])).cast<int>();
      Eigen::Vector3i maxI = Eigen::Vector3d(std::ceil(pos1[0]), std::ceil(pos1[1]), std::ceil(pos1[2])).cast<int>();
      for (int x = minI[0]; x<maxI[0]; x++)
      {
        for (int y = minI[1]; y<maxI[1]; y++)
        {
          for (int z = minI[2]; z<maxI[2]; z++)
          {
            int index = (x-minIndex[0]) + dimensions[0]*(y-minIndex[1]) + dimensions[0]*dimensions[1]*(z-minIndex[2]);
            if (index < 0 || index >= length)
            {
              std::cout << "Error: bad index: " << index << std::endl;
              return;
            }
            // do actual clipping here.... 
            Eigen::Vector3d box_min(((double)x-0.5)*cell_width[0], ((double)y-0.5)*cell_width[1], ((double)z-0.5)*cell_width[2]);
            Cuboid cuboid(box_min, box_min + cell_width);
            Eigen::Vector3d start = starts[i];
            Eigen::Vector3d end = ends[i];

            if (cuboid.clipRay(start, end))
            {
              RGBA col = colours[i];
              if (cells[index].fileName() == "") // first time in this cell, so start writing to a new file
              {
                std::stringstream name;
                name << cloud_name_stub << "_" << x << "_" << y << "_" << z << ".ply"; 
                cells[index].begin(name.str());
              } 
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
        if (cells[i].fileName() == "")
          std::cout << "bad" << std::endl;
        cells[i].writeChunk(chunks[i]);
        chunks[i].clear();
      }
    }       
  };
  if (!Cloud::read(file_name, per_chunk))
    return false;

  for (int i = 0; i<length; i++)
  {
    if (cells[i].fileName() != "")
    {
      cells[i].end();
    }
  }  
  return true;
}

}  // namespace ray
