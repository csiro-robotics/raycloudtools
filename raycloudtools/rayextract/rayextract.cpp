// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/extraction/rayterrain.h"
#include "raylib/extraction/rayforest.h"
#include "raylib/rayutils.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imagewrite.h"

#define _TEST // run multiple times

void usage(bool error=false)
{
  std::cout << "Extract feature into a text file structure" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayextract forest cloud.ply                      - extracts tree locations to file" << std::endl;
  std::cout << "                            --tree_roundness 2   - 1: willow, 0.5: birch, 0.2: pine (length per crown radius)." << std::endl;
  std::cout << "                            --average_height 10  - tree length, if known. -1 to use lowest point as ground instead." << std::endl;
  //  cout << "                             --extrapolate  - estimates tree distribution and adds trees where there is no evidence to the contrary" << endl;
  std::cout << std::endl;
  std::cout << "rayextract terrain cloud.ply             - extract rough terrain undersurface, to mesh." << std::endl;
  std::cout << "                             --width 0.5 - ground width to get average of. Default 0." << std::endl;
  std::cout << std::endl;
  std::cout << "                            --verbose  - extra debug output." << std::endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument file;
  ray::TextArgument forest("forest"), terrain("terrain");
  ray::DoubleArgument tree_roundness(0.01, 3.0);
  ray::DoubleArgument average_height(0.5, 500.0);
  ray::DoubleArgument width(0.0, 10.0);
  ray::OptionalKeyValueArgument roundness_option("tree_roundness", &tree_roundness);
  ray::OptionalKeyValueArgument height_option("average_height", &average_height);
  ray::OptionalKeyValueArgument width_option("width", &width);
  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_forest = ray::parseCommandLine(argc, argv, {&forest, &file}, {&roundness_option, &height_option, &verbose});
  bool extract_terrain = ray::parseCommandLine(argc, argv, {&terrain, &file}, {&width_option, &verbose});
  if (!extract_forest && !extract_terrain)
    usage();
  ray::Cloud cloud;
  if (!cloud.load(file.name()))
    usage(true);

  if (extract_forest)
  {
    ray::Forest forest;
    forest.tree_roundness = tree_roundness.value();
    forest.average_height = average_height.value();
    forest.verbose = verbose.isSet();

#define TEST
#if defined TEST
    const int res = 256;
    const double max_tree_height = (double)res / 8.0;
    
    Eigen::ArrayXXd heightfield = Eigen::ArrayXXd::Constant(res, res, -1e10);
    Eigen::ArrayXXd lowfield = Eigen::ArrayXXd::Constant(res, res, 1e10);
    // now lets give it a base hilly floor
    for (int i = 0; i<res; i++)
    {
      for (int j = 0; j<res; j++)
      {
        double h = std::sin(0.1 * ((double)i + 2.0*double(j)));
        heightfield(i,j) = h + ray::random(-0.25, 0.25);
        lowfield(i,j) = heightfield(i,j);
      }
    }

    int num = 250; // 500
    const double radius_to_height = 0.4; // actual radius to height is 2 * radius_to_height^2
    std::vector<Eigen::Vector3d> ps(num);
    for (int i = 0; i<num; i++)
    {
      double height = max_tree_height * std::pow(2.0, ray::random(-2.0, 0.0)); // i.e. 0.25 to 1 of max_height
      ps[i] = Eigen::Vector3d(ray::random(0.0, (double)res-1.0), ray::random(0.0, (double)res - 1.0), height);
      ps[i][2] += 10.0;
    }
    for (int it = 0; it < 10; it++)
    {
      for (auto &p: ps)
      {
        Eigen::Vector3d shift(0,0,0);
        for (auto &other: ps)
        {
          double radius = radius_to_height * (p[2] + other[2]);
          Eigen::Vector3d dif = p - other;
          dif[2] = 0.0;
          if (dif[0] > res/2.0)
            dif[0] -= res;
          if (dif[1] > res/2.0)
            dif[1] -= res;
          double len_sqr = dif.squaredNorm();
          if (len_sqr > 0.0 && len_sqr < radius*radius)
          {
            double len = std::sqrt(len_sqr);
            shift += (dif/len) * (radius-len);
          }
        }
        p += 0.5*shift;
        p[0] = fmod(p[0] + (double)res, (double)res);
        p[1] = fmod(p[1] + (double)res, (double)res);
      }
    }
    // now make height field
    for (auto &p: ps)
    {
      double radius = radius_to_height * p[2];
      for (int x = (int)(p[0] - radius); x<= (int)(p[0]+radius); x++)
      {
        for (int y = (int)(p[1] - radius); y<= (int)(p[1]+radius); y++)
        {
          double X = (double)x - p[0];
          double Y = (double)y - p[1];
          double mag2 = (double)(X*X + Y*Y);
          if (mag2 <= radius*radius)
          {
            double height = p[2] - 1.0/(0.3*p[2]) * mag2; // (p[2]/2.0)*mag2/(radius*radius);
            int xx = (x + res)%res;
            int yy = (y + res)%res;
            height += ray::random(-1.0, 1.0);
            heightfield(xx, yy) = std::max(heightfield(xx, yy), height);
          }
        }
      }
    }
    // add a noisy function:
    if (1)
    {
      double noise = 5.0;
      const int wid = 80;
      double hs[wid][wid];
      for (int i = 0; i<wid; i++)
      {
        for (int j = 0; j<wid; j++)
          hs[i][j] = ray::random(-noise, noise);
      }
      for (int i = 0; i<res; i++)
      {
        double x = (double)wid * (double)i/(double)res;
        int X = (int)x;
        double blendX = x-(double)X;
        for (int j = 0; j<res; j++)
        {
          double y = (double)wid * (double)j/(double)res;
          int Y = (int)y;
          double blendY = y-(double)Y;
          heightfield(i, j) += hs[X][Y]*(1.0-blendX)*(1.0-blendY) + hs[X][(Y+1)%wid]*(1.0-blendX)*blendY + 
                              hs[(X+1)%wid][Y]*blendX*(1.0-blendY) + hs[(X+1)%wid][(Y+1)%wid]*blendX*blendY; 
          heightfield(i, j) = std::max(heightfield(i, j), -20.0);
        }
      }
    }
    // now render it 
    forest.drawHeightField("highfield.png", heightfield);
    forest.drawHeightField("lowfield.png", lowfield);
    forest.extract(heightfield, lowfield, 1.0);
#else
    forest.extract(cloud);
#endif
  }
  else if (extract_terrain)
  {
    ray::Terrain terrain;
    terrain.extract(cloud, file.name(), width.value(), verbose.isSet());
  }
  else
    usage(true);
}
