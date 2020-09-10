// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raytreegen.h"
#include "raylib/rayforestgen.h"
#include "raylib/rayroomgen.h"
#include "raylib/raybuildinggen.h"
#include "raylib/rayterraingen.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Generates simple example ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycreate room 3 - generates a room using the seed 3. Also:" << std::endl;
  std::cout << "          building" << std::endl;
  std::cout << "          tree" << std::endl;
  std::cout << "          forest" << std::endl;
  std::cout << "          terrain" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::KeyChoice cloud_type({"room", "building", "tree", "forest", "terrain"});
  ray::IntArgument seed(1,1000000);
  if (!ray::parseCommandLine(argc, argv, {&cloud_type, &seed}))
    usage();

  srand(seed.value());

  ray::Cloud cloud;
  const double time_delta = 0.001; // between rays
  std::string type = cloud_type.selectedKey();
  if (type == "room")
  {
    // create room
    ray::RoomGen room_gen;
    room_gen.generate();
    cloud.starts = room_gen.rayStarts();
    cloud.ends = room_gen.rayEnds();
    cloud.times.resize(cloud.starts.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      if (i == cloud.starts.size() / 2)
        time += 0.5;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
    const std::vector<bool> &bounded = room_gen.rayBounded();
    for (int i = 0; i < (int)cloud.colours.size(); i++) 
      cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "building")
  {
    // create building...
    ray::BuildingGen building_gen;
    building_gen.generate();
    cloud.starts = building_gen.rayStarts();
    cloud.ends = building_gen.rayEnds();
    cloud.times.resize(cloud.ends.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
    const std::vector<bool> &bounded = building_gen.rayBounded();
    for (int i = 0; i < (int)cloud.colours.size(); i++) 
      cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "tree" || type == "forest")
  {
    ray::fillBranchAngleLookup();
    double density = 500.0;
    Eigen::Vector3d box_min(-2.0, -2.0, -0.025), box_max(2.0, 2.0, 0.025);
    double time = 0.0;
    if (type == "tree")
    {
      ray::TreeGen tree_gen;
      tree_gen.make(Eigen::Vector3d(0, 0, 0), 0.1, 0.25);
      tree_gen.generateRays(density);
      cloud.starts = tree_gen.rayStarts();
      cloud.ends = tree_gen.rayEnds();
      cloud.times.resize(cloud.ends.size());
      for (size_t i = 0; i < cloud.starts.size(); i++)
      {
        cloud.times[i] = time;
        time += time_delta;
      }
      colourByTime(cloud.times, cloud.colours);
    }
    else if (type == "forest")
    {
      ray::ForestGen forest_gen;
      forest_gen.make(0.25);
      forest_gen.generateRays(density);
      for (auto &tree : forest_gen.trees())
      {  
        const std::vector<Eigen::Vector3d> &ray_starts = tree.rayStarts();
        const std::vector<Eigen::Vector3d> &ray_ends = tree.rayEnds();
        cloud.starts.insert(cloud.starts.end(), ray_starts.begin(), ray_starts.end());
        cloud.ends.insert(cloud.ends.end(), ray_ends.begin(), ray_ends.end());
      }
      cloud.times.resize(cloud.starts.size());
      for (size_t i = 0; i < cloud.starts.size(); i++)
      {
        cloud.times[i] = time;
        time += time_delta;
      }
      box_min *= 2.5;
      box_max *= 2.5;
    }
    int num = int(0.25 * density * (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]));
    for (int i = 0; i < num; i++)
    {
      Eigen::Vector3d pos(ray::random(box_min[0], box_max[0]), ray::random(box_min[1], box_max[1]), ray::random(box_min[2], box_max[2]));
      cloud.ends.push_back(pos);
      cloud.starts.push_back(pos + Eigen::Vector3d(ray::random(-0.1, 0.1), ray::random(-0.1, 0.1), ray::random(0.2, 0.5)));
      cloud.times.push_back(time);
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else if (type == "terrain")
  {
    ray::TerrainGen terrain;
    terrain.generate();
    cloud.starts = terrain.rayStarts();
    cloud.ends = terrain.rayEnds();
    cloud.times.resize(cloud.starts.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else
    usage();
  cloud.save(type + ".ply");

  return true;
}