// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Split a ray cloud relative to the supplied triangle mesh, generating two cropped ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysplit raycloud pos 10,0,0             - splits along x axis" << std::endl;
  std::cout << "                  time 10000             - splits at given acquisition time" << std::endl;
  std::cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << std::endl;
  std::cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << std::endl;
  std::cout << "                  meshfile distance 0.2  - splits raycloud at 0.2m from the meshfile surface" << std::endl;
  std::cout << "                  startpos 1,2,3         - splits based on start position, around plane 1,2,3" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays"
       << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  speed 1.0              - splits out rays when sensor moving above the given speed" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::Vector3dArgument pos, colour(0.0, 1.0), startpos, raydir(-1.0, 1.0);
  ray::DoubleArgument time, alpha(0.0,1.0), range(0.0,1000.0), speed(0.0,1000.0);
  ray::KeyValueChoice choice({"pos", "time", "colour", "alpha", "startpos", "raydir", "range", "speed"}, 
                             {&pos,  &time,  &colour,  &alpha,  &startpos,  &raydir,  &range,  &speed});
  ray::FileArgument mesh_file;
  ray::TextArgument text("distance");
  ray::DoubleArgument mesh_offset;
  bool standard_format = ray::parseCommandLine(argc, argv, {&cloud_file, &choice});
  bool mesh_split = ray::parseCommandLine(argc, argv, {&cloud_file, &mesh_file, &text, &mesh_offset});
  if (!standard_format && !mesh_split)
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();
  ray::Cloud inside, outside;
  if (mesh_split)
  {
    ray::Mesh mesh;
    ray::readPlyMesh(mesh_file.name(), mesh);
    mesh.splitCloud(cloud, mesh_offset.value(), inside, outside);
  }
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      cloud.split(inside, outside, [&](int i) -> bool { return cloud.times[i] > time.value(); });
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      cloud.split(inside, outside, [&](int i) { return cloud.colours[i].alpha > c; });
    }
    else if (parameter == "pos")
    {
      Eigen::Vector3d vec = pos.value() / pos.value().squaredNorm();
      cloud.split(inside, outside, [&](int i) { return cloud.ends[i].dot(vec) > 1.0; });
    }
    else if (parameter == "startpos")
    {
      Eigen::Vector3d vec = startpos.value() / startpos.value().squaredNorm();
      cloud.split(inside, outside, [&](int i) { return cloud.starts[i].dot(vec) > 0.0; });
    }
    else if (parameter == "raydir")
    {
      Eigen::Vector3d vec = raydir.value() / raydir.value().squaredNorm();
      cloud.split(inside, outside, [&](int i) {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 0.0;
      });
    }
    else if (parameter == "colour")
    {
      Eigen::Vector3d vec = colour.value() / colour.value().squaredNorm();
      cloud.split(inside, outside, [&](int i) {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                     (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 0.0;
      });
    }
    else if (parameter == "range")
    {
      cloud.split(inside, outside, [&](int i) { return (cloud.starts[i] - cloud.ends[i]).norm() > range.value(); });
    }
    else if (parameter == "speed")
    {
      cloud.split(inside, outside, [&](int i) {
        if (i == 0)
          return false;
        return (cloud.starts[i] - cloud.starts[i - 1]).norm() / (cloud.times[i] - cloud.times[i - 1]) > speed.value();
      });
    }
  }

  std::string file_stub = cloud_file.nameStub();
  inside.save(file_stub + "_inside.ply");
  outside.save(file_stub + "_outside.ply");
  return true;
}
