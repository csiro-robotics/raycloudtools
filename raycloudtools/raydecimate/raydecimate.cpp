// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 1)
{
  std::cout << "Decimate a ray cloud spatially or temporally" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << std::endl;
  std::cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::IntArgument num_rays(1, 100);
  ray::DoubleArgument vox_width(0.01, 100.0);
  ray::ValueKeyChoice quantity({ &vox_width, &num_rays }, { "cm", "rays" });
  if (!ray::parseCommandLine(argc, argv, { &cloud_file, &quantity }))
    usage();
  const bool spatial_decimation = quantity.selectedKey() == "cm";

  ray::CloudWriter writer;
  if (!writer.begin(cloud_file.nameStub() + "_decimated.ply"))
    usage();

  // By maintaining these buffers below, we avoid almost all memory fragmentation
  ray::Cloud chunk;
  std::vector<int64_t> subsample;
  // voxel set is global, however its size is proportional to the decimated cloud size,
  // so we expect it to fit within RAM limits
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;

  auto decimate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    if (spatial_decimation)
    {
      double width = 0.01 * vox_width.value();
      subsample.clear();
      voxelSubsample(ends, width, subsample, voxel_set);
      chunk.resize(subsample.size());
      for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
      {
        int64_t id = subsample[i];
        chunk.starts[i] = starts[id];
        chunk.ends[i] = ends[id];
        chunk.colours[i] = colours[id];
        chunk.times[i] = times[id];
      }
    }
    else
    {
      size_t decimation = (size_t)num_rays.value();
      size_t count = (ends.size() + decimation - 1) / decimation;
      chunk.resize(count);
      for (size_t i = 0, c = 0; i < ends.size(); i += decimation, c++)
      {
        chunk.starts[c] = starts[i];
        chunk.ends[c] = ends[i];
        chunk.times[c] = times[i];
        chunk.colours[c] = colours[i];
      }
    }
    writer.writeChunk(chunk);
  };

  if (!ray::Cloud::read(cloud_file.name(), decimate))
    usage();
  writer.end();

  return 0;
}
