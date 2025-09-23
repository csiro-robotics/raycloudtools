// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"

#include <nabo/nabo.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imageread.h"
#include "raylib/imagewrite.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Smooth a ray cloud. Nearby off-surface points are moved onto the nearest surface." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysmooth raycloud" << std::endl;
  // clang-format on
  exit(exit_code);
}

int raySmooth(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  if (!ray::parseCommandLine(argc, argv, { &cloud_file }))
    usage();

  #define INPAINT
  #if defined INPAINT
  int width, height, bytes;
  unsigned char *data = stbi_load(cloud_file.name().c_str(), &width, &height, &bytes, 0);
  std::vector<ray::RGBA> pixel_colours(width * height, ray::RGBA(0,0,0,0));

  for (int x = 0; x<width; x++)
  {
    if (!(x%100))
      std::cout << "computed " << x << " / " << width << std::endl;
    for (int y = 0; y<height; y++)
    {
      for (int r = 0; r<3; r++)
      {
        int J = (x + width*y)*4;
        if (data[J+3])
        {
          pixel_colours[x + width*y] = *(ray::RGBA *)&data[J];
          break;
        }
        Eigen::Vector3d colour(0,0,0);
        double count = 0.0;
        for (int i = -r; i<=r; i++)
        {
          int xs[4] = {x+i, x+i, x-r, x+r};
          int ys[4] = {y-r, y+r, y+i, y+i};
          for (int id = 0; id<4; id++)
          {
            if (xs[id] < 0 || xs[id] >= width || ys[id] < 0 || ys[id] >= height)
              continue;
            int I = (xs[id] + width*ys[id])*4;
            if (data[I + 3])
            {
              ray::RGBA col = *(ray::RGBA *)&data[I];
              colour += Eigen::Vector3d(col.red, col.green, col.blue);
              count++;
            }
          }
        }
        if (count > 0)
        {
          colour /= count;
          pixel_colours[x + width*y] = ray::RGBA((uint8_t)colour[0], (uint8_t)colour[1], (uint8_t)colour[2], 255);
          break;
        }
      }
    }
  }

  stbi_image_free(data);
 // stbi_flip_vertically_on_write(1);
  stbi_write_png((cloud_file.nameStub() + "_inpainted.png").c_str(), width, height, 4, (void *)&pixel_colours[0], 4 * width);

  return 1;
  #endif
  #define TRANSFORM_TRAJ
  #if defined TRANSFORM_TRAJ
  std::ifstream ifs(cloud_file.name(), std::ios::in);
  std::string line;
  if (!ifs)
  {
    std::cerr << "Failed to open trajectory file" << std::endl;
    return false;
  }
  std::vector<Eigen::Vector3d> points;
  std::vector<double> times;
  std::vector<Eigen::Quaterniond> quats;
  while (ifs.good())
  {
    if (ifs.fail())
    {
      std::cerr << "Invalid stream when loading trajectory file" << std::endl;
      return false;
    }
    getline(ifs, line);
    if (line.length() == 0 || line[0] == '%')
      continue;

    std::istringstream iss(line);
    double time;
    Eigen::Vector3d point;
    Eigen::Quaterniond quat;
    iss >> time >> point[0] >> point[1] >> point[2] >> quat.w() >> quat.x() >> quat.y() >> quat.z();
    if (iss.fail())
    {
      std::cerr << "Invalid fields" << std::endl;
      return false;
    }
    times.push_back(time);

   // #define REALIGN_HANDHELD
    #if defined REALIGN_HANDHELD
      point[0] += 8.0;
      Eigen::AngleAxisd a(4.5*(ray::kPi/180.0), Eigen::Vector3d(0,0,1));
      quat = Eigen::Quaterniond(a) * quat;
    #endif

    points.push_back(point);
    quats.push_back(quat);
  }  

  for (int j = 0; j<4; j++)
  {
    std::ofstream ofs((cloud_file.nameStub() + std::to_string(j) + ".txt").c_str(), std::ios::out);
    ofs << "\%time x y z qw qx qy qz userfields\n";
    for (int i = 0; i<(int)points.size(); i++)
    {
      Eigen::AngleAxisd aa(0.5*ray::kPi * (double)j, Eigen::Vector3d(0,0,1));
      Eigen::Quaterniond rot(aa);
      Eigen::Quaterniond quat = rot * quats[i];
      ofs << std::setprecision(16) << times[i] << " " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << quat.w() << " " << quat.x() << " " << quat.y() << " " << quat.z() << "\n";
    }
    ofs.close();
  }
  return 1;
  #endif
  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();

  // Method:
  // 1. generate normals and neighbour indices
  // 2. pull point along normal direction so as to match neighbours, weighted by normal similarity

  const int num_neighbours = 16;
  std::vector<Eigen::Vector3d> normals;
  Eigen::MatrixXi neighbour_indices;
  cloud.getSurfels(num_neighbours, nullptr, &normals, nullptr, nullptr, &neighbour_indices);

  std::vector<Eigen::Vector3d> centroids(cloud.ends.size());
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    double total_weight = 0.2;  // more averaging if it uses less of the central position, but 0 risks a divide by 0
    Eigen::Vector3d weighted_sum = cloud.ends[i] * total_weight;
    for (int j = 0; j < num_neighbours && neighbour_indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
    {
      int k = neighbour_indices(j, i);
      double weight = std::max(0.0, 1.0 - (normals[k] - normals[i]).squaredNorm());
      weighted_sum += cloud.ends[k] * weight;
      total_weight += weight;
    }
    centroids[i] = weighted_sum / total_weight;
  }
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    cloud.ends[i] += normals[i] * (centroids[i] - cloud.ends[i]).dot(normals[i]);
  }

  cloud.save(cloud_file.nameStub() + "_smooth.ply");

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(raySmooth, argc, argv);
}