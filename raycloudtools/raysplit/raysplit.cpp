// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"

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
  if (argc != 4 && argc != 5)
    usage();

  std::string file = argv[1];
  ray::Cloud cloud;
  cloud.load(file);

  ray::Cloud inside, outside;
  if (argc == 5)
  {
    std::string mesh_file = argv[2];
    ray::Mesh mesh;
    ray::readPlyMesh(mesh_file, mesh);

    double offset = std::stod(argv[4]);
    mesh.splitCloud(cloud, offset, inside, outside);
  }
  else
  {
    std::string parameter = std::string(argv[2]);
    if (parameter == "time")
    {
      double val = std::stod(argv[3]);
      cloud.split(inside, outside, [&](int i) -> bool { return cloud.times[i] > val; });
    }
    else if (parameter == "alpha")
    {
      double val = std::stod(argv[3]);
      if (!(val >= 0.0 && val <= 1.0))
        usage();
      uint8_t c = uint8_t(255.0 * val);
      cloud.split(inside, outside, [&](int i) { return cloud.colours[i].alpha > c; });
    }
    else if (parameter == "pos")
    {
      std::stringstream ss(argv[3]);
      Eigen::Vector3d vec;
      ss >> vec[0];
      ss.ignore(1);
      ss >> vec[1];
      ss.ignore(1);
      ss >> vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, [&](int i) { return cloud.ends[i].dot(vec) > 1.0; });
    }
    else if (parameter == "startpos")
    {
      std::stringstream ss(argv[3]);
      Eigen::Vector3d vec;
      ss >> vec[0];
      ss.ignore(1);
      ss >> vec[1];
      ss.ignore(1);
      ss >> vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, [&](int i) { return cloud.starts[i].dot(vec) > 0.0; });
    }
    else if (parameter == "raydir")
    {
      std::stringstream ss(argv[3]);
      Eigen::Vector3d vec;
      ss >> vec[0];
      ss.ignore(1);
      ss >> vec[1];
      ss.ignore(1);
      ss >> vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, [&](int i) {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 0.0;
      });
    }
    else if (parameter == "colour")
    {
      std::stringstream ss(argv[3]);
      Eigen::Vector3d vec;
      ss >> vec[0];
      ss.ignore(1);
      ss >> vec[1];
      ss.ignore(1);
      ss >> vec[2];
      vec /= vec.squaredNorm();

      cloud.split(inside, outside, [&](int i) {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                     (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 0.0;
      });
    }
    else if (parameter == "range")
    {
      double val = std::stod(argv[3]);
      cloud.split(inside, outside, [&](int i) { return (cloud.starts[i] - cloud.ends[i]).norm() > val; });
    }
    else if (parameter == "speed")
    {
      double val = std::stod(argv[3]);
      cloud.split(inside, outside, [&](int i) {
        if (i == 0)
          return false;
        return (cloud.starts[i] - cloud.starts[i - 1]).norm() / (cloud.times[i] - cloud.times[i - 1]) > val;
      });
    }
  }

  std::string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);

  inside.save(file_stub + "_inside.ply");
  outside.save(file_stub + "_outside.ply");
  return true;
}
