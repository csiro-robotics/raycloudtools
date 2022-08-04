// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayroomgen.h"
namespace ray
{
// A room with a door, window, table and cupboard
void RoomGen::generate()
{
  double point_density = 750.0;
  double room_width = random(3.0, 6.0);
  double room_length = random(3.0, 6.0);
  double room_height = random(2.75, 3.0);

  Eigen::Vector3d floor_centre(0, 0,
                               -room_height * 0.5);  // random(-10.0, 10.0), random(-10.0, 10.0), random(-10.0, 10.0));
  double room_yaw = random(0.0, 2.0 * kPi);

  std::vector<Cuboid> negatives;
  Cuboid room(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(room_width, room_length, room_height));
  negatives.push_back(room);

  double door_width = 0.7;
  double door_start = random(0.0, room_width - door_width);
  double door_height = random(2.2, 2.5);
  Eigen::Vector3d door_pos = Eigen::Vector3d(door_start, -0.2, 0.0);
  Cuboid door(door_pos, door_pos + Eigen::Vector3d(door_width, 0.3, door_height));
  negatives.push_back(door);

  Cuboid outside_door(Eigen::Vector3d(-20.0, -20.0, 0.0), Eigen::Vector3d(20.0, -0.15, 20.0));
  negatives.push_back(outside_door);

  double window_width = random(0.7, 1.5);
  double window_height = random(0.9, 1.4);
  double window_start = random(0.4, room_length - window_width - 0.4);
  Eigen::Vector3d window_pos = Eigen::Vector3d(-0.2, window_start, 1.0);
  Cuboid window(window_pos, window_pos + Eigen::Vector3d(0.3, window_width, window_height));
  negatives.push_back(window);

  Cuboid outside_window(Eigen::Vector3d(-20.0, -20.0, 0.0), Eigen::Vector3d(-0.15, 20.0, 20.0));
  negatives.push_back(outside_window);

  std::vector<Cuboid> positives;
  double table_width = random(0.5, 1.5);
  double table_length = random(0.5, 1.5);
  double table_height = random(0.5, 1.2);
  Eigen::Vector3d table_pos(random(0.0, room_width - table_width - 1.0), random(0.0, room_length - table_length - 1.0),
                            table_height);
  Cuboid table_top(table_pos, table_pos + Eigen::Vector3d(table_width, table_length, 0.05));
  positives.push_back(table_top);
  for (int x = 0; x < 2; x++)
  {
    for (int y = 0; y < 2; y++)
    {
      Eigen::Vector3d pos(table_pos[0] + 0.05 + (table_width - 0.15) * (double)x,
                          table_pos[1] + 0.05 + (table_length - 0.15) * (double)y, 0.0);
      Cuboid leg(pos, pos + Eigen::Vector3d(0.05, 0.05, table_height));
      positives.push_back(leg);
    }
  }

  double cupboard_width = random(1.0, 2.9);
  double cupboard_start = random(0.0, room_length - cupboard_width);
  double cupboard_depth = random(0.2, 0.8);
  double cupboard_height = random(1.0, 2.5);
  Eigen::Vector3d cupboard_pos(room_width - cupboard_depth, cupboard_start, 0.0);
  Cuboid cupboard(cupboard_pos, cupboard_pos + Eigen::Vector3d(cupboard_depth, cupboard_width, cupboard_height));
  positives.push_back(cupboard);


  // OK now we ray trace from some random location onto the set of cuboids...
  Eigen::Vector3d start(random(1.4, room_width - 1.4), random(1.4, room_length - 1.4), random(1.3, 2.0));
  Eigen::Vector3d s = start - Eigen::Vector3d(room_width, room_length, 0.0) / 2.0;
  Eigen::Vector3d ray_start =
    Eigen::Vector3d(s[0] * cos(room_yaw) + s[1] * sin(room_yaw), -s[0] * sin(room_yaw) + s[1] * cos(room_yaw), s[2]) +
    floor_centre;
  size_t num_rays = (size_t)(point_density * (room_width * room_length) * 2.0 + room_width * room_height * 2.0 +
                             room_length * room_height * 2.0);
  for (size_t i = 0; i < num_rays; i++)
  {
    Eigen::Vector3d dir(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0));
    dir.normalize();
    const double max_range = 20.0;
    std::vector<Eigen::Vector3d> hits;
    for (int i = 0; i < (int)negatives.size(); i++)
    {
      double new_range = max_range;
      if (negatives[i].intersectsRay(start, dir, new_range, false))
        hits.push_back(start + dir * (new_range + 1e-6));
    }
    double range = max_range;
    for (auto &hit : hits)
    {
      bool intersected = false;
      for (auto &cuboid : negatives)
      {
        if (cuboid.intersects(hit))
        {
          intersected = true;
          break;
        }
      }
      if (!intersected)
        range = std::min(range, (hit - start).norm());
    }
    if (i > num_rays / 2)
    {
      for (auto &cuboid : positives) cuboid.intersectsRay(start, dir, range, true);
    }

    const double range_noise = 0.03;
    Eigen::Vector3d end = start + (range + random(-range_noise, range_noise)) * dir;
    Eigen::Vector3d s = end - Eigen::Vector3d(room_width, room_length, 0.0) / 2.0;
    Eigen::Vector3d ray_end =
      Eigen::Vector3d(s[0] * cos(room_yaw) + s[1] * sin(room_yaw), -s[0] * sin(room_yaw) + s[1] * cos(room_yaw), s[2]) +
      floor_centre;
    ray_starts_.push_back(ray_start);
    ray_ends_.push_back(ray_end);
    ray_bounded_.push_back(range != max_range);
  }
}
}  // namespace ray