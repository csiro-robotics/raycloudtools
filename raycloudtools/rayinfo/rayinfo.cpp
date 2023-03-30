// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/rayply.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <chrono>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Get basic information on the ray cloud, such as its bounds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayinfo raycloud.ply" << std::endl;
  // clang-format on
  exit(exit_code);
}

std::string getTime(double unix_time)
{
  time_t unix_timestamp = static_cast<time_t>(unix_time);
  char time_buf[80];
  struct tm ts;
  ts = *localtime(&unix_timestamp);
  strftime(time_buf, sizeof(time_buf), "%a %d/%m/%Y %H:%M:%S %Z", &ts);

  return std::string(time_buf);
}


int main(int argc, char *argv[])
{
  ray::FileArgument cloud;
  if (!ray::parseCommandLine(argc, argv, { &cloud }, {}))
  {
    usage();
  }

  ray::Cloud::Info info;
  ray::Cloud::getInfo(cloud.name(), info);

  // path length
  // num jumps
  int out_of_order = 0;
  int time_jumps = 0;
  int space_jumps = 0;
  int num_jumps = 0;
  double last_time = std::numeric_limits<double>::lowest();
  double path_length = 0.0, path_period = 0.0;
  Eigen::Vector3d last_start(0,0,0);
  double min_ray_length = std::numeric_limits<double>::max();
  double max_ray_length = std::numeric_limits<double>::lowest();
  double max_bounded_ray_length = std::numeric_limits<double>::lowest();
  ray::RGBA min_col(255, 255, 255, 255), max_col(0, 0, 0, 0);
  auto get_info = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (last_time != std::numeric_limits<double>::lowest())
      {
        if (times[i] < last_time)
        {
          out_of_order++;
        }
        double distance = (starts[i] - last_start).norm();
        double time_delta = std::abs(times[i] - last_time);
        bool time_jump = time_delta > 1.00001;
        bool space_jump = distance > 1.00001;
        if (time_jump || space_jump)
        {
          time_jumps += time_jump ? 1 : 0;
          space_jumps += space_jump ? 1 : 0;
          num_jumps++;
        }
        else
        {
          path_length += distance;
          path_period += time_delta;
        }
      }
      last_start = starts[i];
      last_time = times[i];

      double ray_length = (ends[i] - starts[i]).norm();
      min_ray_length = std::min(min_ray_length, ray_length);
      max_ray_length = std::max(max_ray_length, ray_length);
      if (colours[i].alpha > 0)
      {
        max_bounded_ray_length = std::max(max_bounded_ray_length, ray_length);

        min_col.red = std::min(min_col.red, colours[i].red);
        min_col.green = std::min(min_col.green, colours[i].green);
        min_col.blue = std::min(min_col.blue, colours[i].blue);
        min_col.alpha = std::min(min_col.alpha, colours[i].alpha);
        max_col.red = std::max(max_col.red, colours[i].red);
        max_col.green = std::max(max_col.green, colours[i].green);
        max_col.blue = std::max(max_col.blue, colours[i].blue);
        max_col.alpha = std::max(max_col.alpha, colours[i].alpha);
      }
      else
      {
        min_col.alpha = 0;
      }
    }
  };
  if (!ray::readPly(cloud.name(), true, get_info, 0))
  {
    usage();
  }

  std::cout << std::endl;
  std::cout << "Ray cloud information for " << cloud.name() << ":" << std::endl;
  std::cout << std::endl;
  std::cout << "  number of rays: \t" << info.num_bounded + info.num_unbounded << " of which " << info.num_bounded << " have end points." << std::endl;
  int min_ms = (int)info.min_time/1e6;
  info.min_time -= 1e6 * (double)min_ms;
  int max_ms = (int)info.max_time/1e6;
  info.max_time -= 1e6 * (double)max_ms;
  std::cout << "  date from: \t\t" << getTime(info.min_time) << "\t(" << min_ms << " Ms \t" << info.min_time << " s)" << std::endl;
  std::cout << "         to: \t\t" << getTime(info.max_time) << "\t(" << max_ms << " Ms \t" << info.max_time << " s)" << std::endl;
  std::cout << "  full bounds: \t\t" << info.rays_bound.min_bound_.transpose() << " to " << info.rays_bound.max_bound_.transpose() << std::endl;
  std::cout << "  bounds of end points:\t" << info.ends_bound.min_bound_.transpose() << " to " << info.ends_bound.max_bound_.transpose() << std::endl;
  std::cout << "  first location: \t" << info.start_pos.transpose() << ", last location: " << info.end_pos.transpose() << std::endl;
  std::cout << std::endl;
  if (out_of_order > (info.num_bounded + info.num_unbounded)/16)
  {
    std::cout << "  times are unordered." << std::endl;
  }
  else if (num_jumps > (info.num_bounded + info.num_unbounded)/16)
  {
    std::cout << "  ray starts are discontinuous." << std::endl;
  }
  else
  {
    std::cout << "  contiguous blocks: \t" << num_jumps+1;
    if (num_jumps > 0)
    {
      std::cout << " (discontinuities: " << time_jumps << " > 1s, " << space_jumps << " > 1m)";
    }
    std::cout << std::endl;
    std::cout << "  path length: \t\t" << path_length << " m and period: " << path_period << " s" << std::endl;
  }
  std::cout << "  ray length: \t\t" << min_ray_length << " to " << max_ray_length << " m, max end point ray length: " << max_bounded_ray_length << " m" << std::endl;
  std::cout << "  colour range (RGBA): \t" << (int)min_col.red << "," << (int)min_col.green << "," << (int)min_col.blue << "," << (int)min_col.alpha << " to " << 
                                            (int)max_col.red << "," << (int)max_col.green << "," << (int)max_col.blue << "," << (int)max_col.alpha << std::endl;
 
  return 0;
}
