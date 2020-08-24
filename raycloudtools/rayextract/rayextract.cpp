// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/extraction/rayterrain.h"
#include "raylib/extraction/rayforest.h"
#include "raylib/rayutils.h"

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
  std::cout << "rayextract forest cloud.ply --ground_mesh mesh.ply - extracts tree locations to file, using a known ground mesh" << std::endl;
  std::cout << "                            --tree_roundness 2   - 1: willow, 0.5: birch, 0.2: pine (length per crown radius)." << std::endl;
  std::cout << "                            --average_height 10       - use when heights are uniform, shapes can vary." << std::endl;
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
  if (argc < 3 || argc > 12)
    usage(true);
  std::string type(argv[1]);
  ray::Cloud cloud;
  std::string file(argv[2]);
  if (!cloud.load(file))
    usage(true);
  if (file.substr(file.length() - 4) == ".ply")
    file = file.substr(0, file.length() - 4);

  if (type == "forest")
  {
    double roundness = 0.0;
    double height = 0.0;
    bool verbose = false;
    for (int i = 3; i<argc; i+=2)
    {
      std::string str(argv[i]);
      if (str == "--tree_roundness" || str == "-t")
        roundness = std::stod(argv[i+1]);
      else if (str == "--average_height" || str == "-a")
        height = std::stod(argv[i+1]);
      else if (str == "--verbose" || str == "-v")
        verbose = true;
    }
    ray::Forest forest;
    forest.tree_roundness = roundness;
    forest.average_height = height;
    forest.verbose = verbose;

#define TEST
#if defined TEST
    const int res = 256;
    const double max_tree_height = (double)res / 8.0;
    
    Eigen::ArrayXXd heightfield(res, res);
    memset(&heightfield(0,0), 0, res*res*sizeof(double));

    int num = 500;
    const double radius_to_height = 0.4;
    std::vector<Eigen::Vector3d> ps(num);
    for (int i = 0; i<num; i++)
    {
      double height = max_tree_height * std::pow(2.0, ray::random(-2.0, 0.0)); // i.e. 0.25 to 1 of max_height
      ps[i] = Eigen::Vector3d(ray::random(0.0, (double)res-1.0), ray::random(0.0, (double)res - 1.0), height);
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
            double height = p[2] - (p[2]/2.0)*mag2/(radius*radius);
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
      const int wid = 80;
      double hs[wid][wid];
      for (int i = 0; i<wid; i++)
      {
        for (int j = 0; j<wid; j++)
          hs[i][j] = ray::random(-5.0, 5.0);
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
          heightfield(i, j) += hs[X][Y]*(1.0-blendX)*(1.0-blendY) + hs[X][Y+1]*(1.0-blendX)*blendY + 
                              hs[X+1][Y]*blendX*(1.0-blendY) + hs[X+1][Y+1]*blendX*blendY; 
          heightfield(i, j) = std::max(heightfield(i, j), 0.0);
        }
      }
    }
    // now render it 
    forest.drawHeightField("testheight.png", heightfield);
    forest.extract(heightfield, 1.0);
#else
    forest.extract(cloud);
#endif
  }
  else if (type == "terrain")
  {
    bool verbose = false;
    double width = 0.0;
    for (int i = 3; i<argc; i++)
    {
      if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v")
        verbose = true;
      if (std::string(argv[i]) == "--width" || std::string(argv[i]) == "-m")
      {
        if (argc == i+1)
          usage(true);
        width = std::stod(argv[i+1]);
      }
    }
    ray::Terrain terrain;
    terrain.extract(cloud, file, width, verbose);
  }
  else
    usage(true);
}
