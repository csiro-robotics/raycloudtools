// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Decimate a ray cloud spatially or temporally" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydecimate raycloud 3 cm   - reduces to one end point every 3 cm" << std::endl;
  std::cout << "raydecimate raycloud 4 rays - reduces to every fourth ray" << std::endl;
  std::cout << "raydecimate raycloud 0.1 seconds - reduces to one ray every 0.1 seconds" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 4)
    usage();

  std::string file = argv[1];
  ray::Cloud cloud;
  cloud.load(file);
  std::vector<Eigen::Vector3d> starts, ends;
  std::vector<ray::RGBA> colours;
  std::vector<double> times;
  std::string type = argv[3];
  if (type == "cm")
  {
    cloud.decimate(0.01 * std::stod(argv[2]));
  }
  else if (type == "rays")
  {
    int decimation = std::stoi(argv[2]);
    for (int i = 0; i < (int)cloud.ends.size(); i += decimation)
    {
      starts.push_back(cloud.starts[i]);
      ends.push_back(cloud.ends[i]);
      times.push_back(cloud.times[i]);
      colours.push_back(cloud.colours[i]);
    }
    cloud.starts = starts;
    cloud.ends = ends;
    cloud.times = times;
    cloud.colours = colours;
  }
  else if (type == "seconds" || type == "s")
  {
    double delta = std::stod(argv[2]);
    double last_time = cloud.times[0];
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (cloud.times[i] >= last_time + delta)
      {
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        colours.push_back(cloud.colours[i]);
        while (cloud.times[i] >= last_time + delta)  // in case delta is tiny
          last_time += delta;
      }
    }
    cloud.starts = starts;
    cloud.ends = ends;
    cloud.times = times;
    cloud.colours = colours;
  }
  else
    usage(false);
  std::string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);
  cloud.save(file_stub + "_decimated.ply");
  return true;
}
