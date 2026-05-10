// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raybuildinggen.h"
#include "raylib/raycloud.h"
#include "raylib/rayforestgen.h"
#include "raylib/rayparse.h"
#include "raylib/rayroomgen.h"
#include "raylib/rayterraingen.h"
#include "raylib/raytreegen.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Generates simple example ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycreate room 3 - generates a room using the seed 3. Also:" << std::endl;
  std::cout << "          building" << std::endl;
  std::cout << "          tree" << std::endl;
  std::cout << "          forest" << std::endl;
  std::cout << "          terrain" << std::endl;
  std::cout << "          field" << std::endl;
  std::cout << std::endl;
  std::cout << "          forest trees.txt - generate from a comma-separated list of x,y,z,radius trees" << std::endl;
  std::cout << "          terrain mesh.ply      - generate from a ground mesh" << std::endl;
  // clang-format on
  exit(exit_code);
}

double random(double x, double y){ return ray::random(x, y); }

int rayCreate(int argc, char *argv[])
{
  ray::KeyChoice cloud_type({ "room", "building", "tree", "forest", "terrain", "field" });
  ray::IntArgument seed(1, 1000000);
  ray::FileArgument input_file;
  bool from_seed = ray::parseCommandLine(argc, argv, { &cloud_type, &seed });
  bool from_file = ray::parseCommandLine(argc, argv, { &cloud_type, &input_file });
  if (!from_seed && !from_file)
    usage();

  if (from_seed)
  {
    ray::srand(seed.value());
  }

  ray::Cloud cloud;
  const double time_delta = 0.001;  // between rays
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
    for (int i = 0; i < (int)cloud.colours.size(); i++) cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "field") // like a field of wheat. The rays follow a Beer/Lambert law in their depth
  {
    double half_width = random(5.0, 20.0);
    double crop_height = random(0.2, 1.5);
    double height_gradient = 0.0;
    double density_gradient = 0.0;
    double crop_density = random(0.04, 0.5);
    double ray_density = random(400.0, 1600.0); // rays per square metre
    std::cout << "Creating field of width " << 2.0*half_width << " m, height: " << crop_height << " m, crop density " << crop_density << " and " << ray_density << " vertical rays per m^2" << std::endl;
    int num_rays = (int)(ray_density * 4.0 * half_width*half_width);
    if (random(0.0,1.0) < 0.5 && seed.value() != 1)
    {
      height_gradient = random(-crop_height/half_width, crop_height/half_width);
      std::cout << "height gradient in y: " << height_gradient << std::endl;
    }
    if (random(0.0, 1.0) < 0.5)
    {
      density_gradient = random(-crop_density/half_width, crop_density/half_width);  
      std::cout << "density gradient in x: " << density_gradient << std::endl;
    }

    for (int i = 0; i<num_rays; i++)
    {
      Eigen::Vector3d pos(random(-half_width, half_width), random(-half_width, half_width), 0.0);
      double lambda = crop_density + pos[0]*density_gradient;
      double height = crop_height + pos[1]*height_gradient;
      Eigen::Vector3d start(pos[0], pos[1], 2.0*height);

      double prob_ground = 1.0 - std::exp(-lambda * height);
      double rnd = random(0.0, 1.0);
      if (rnd < prob_ground)
      {
        double h = std::log(1.0 - rnd)/-lambda;
        pos[2] = height - h;
        if (pos[2] < -0.01)
          std::cout << "error: " << pos[2] << ", rnd: " << rnd << ", h: " << h << std::endl;
      }
      cloud.addRay(start, pos, 0.001*(double)i, ray::RGBA(255,255,255,255));
    }
    colourByTime(cloud.times, cloud.colours);
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
    for (int i = 0; i < (int)cloud.colours.size(); i++) cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "tree" || type == "forest")
  {
    const double density = 500.0;                   // density of points on the branches of the trees
    const double tree_ground_extent = 2.0;          // extent (from centre) is the half-width
    const double ground_noise_extent = 0.025;       // vertical noise in the ground
    const double forest_ground_multiplier = 5.0;    // how much larger the forest ground is
    const double ground_ray_deviation = 0.1;        // lateral deviation in start of ray, from end point
    const double ground_ray_vertical_height = 1.5;  // height above ground for ray start

    ray::fillBranchAngleLookup();
    Eigen::Vector3d box_min(-tree_ground_extent, -tree_ground_extent, -ground_noise_extent);
    Eigen::Vector3d box_max(tree_ground_extent, tree_ground_extent, ground_noise_extent);
    double time = 0.0;
    if (type == "tree")  // create a single tree
    {
      ray::TreeGen tree_gen;
      ray::TreeParams params;
      params.random_factor = 0.25;
      tree_gen.segments().resize(1);
      tree_gen.segments()[0].tip = Eigen::Vector3d(0, 0, 0);
      tree_gen.segments()[0].radius = 0.1;
      tree_gen.make(params);
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
    else if (type == "forest")  // create multiple trees on a plane
    {
      ray::ForestParams params;
      params.random_factor = 0.25;
      ray::ForestGen forest_gen;
      if (from_file)  // load a forest from an _trees.txt file
      {
        if (!forest_gen.makeFromFile(input_file.name(), params))
        {
          usage();
        }
      }
      else
      {
        forest_gen.make(params);
      }
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
      box_min *= forest_ground_multiplier;  // for a forest, we need a larger ground
      box_max *= forest_ground_multiplier;
    }
    if (!from_file)
    {
      int num = int(0.25 * density * (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]));
      for (int i = 0; i < num; i++)
      {
        Eigen::Vector3d pos(ray::random(box_min[0], box_max[0]), ray::random(box_min[1], box_max[1]),
                            ray::random(box_min[2], box_max[2]));
        cloud.ends.push_back(pos);
        cloud.starts.push_back(pos + Eigen::Vector3d(ray::random(-ground_ray_deviation, ground_ray_deviation),
                                                     ray::random(-ground_ray_deviation, ground_ray_deviation),
                                                     ground_ray_vertical_height));
        cloud.times.push_back(time);
        time += time_delta;
      }
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else if (type == "terrain")
  {
    ray::TerrainGen terrain;
    if (from_file)  // generate ray cloud terrain from a .ply mesh file
    {
      terrain.generateFromFile(input_file.name());
    }
    else
    {
      terrain.generate();
    }
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

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayCreate, argc, argv);
}