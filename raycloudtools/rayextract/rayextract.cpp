// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/extraction/raytrunks.h"
#include "raylib/extraction/rayterrain.h"
#include "raylib/extraction/rayforest.h"
#include "raylib/extraction/raytrees.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayparse.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/extraction/rayclusters.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(bool error=false)
{
  std::cout << "Extract feature into a text file structure" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayextract trunks cloud.ply                 - estimate tree trunks and save to text file" << std::endl;
  std::cout << "rayextract trees cloud.ply cloud_trunks.txt - estimate trees using trunks as seeds, and save to text file" << std::endl;
  std::cout << "rayextract forest cloud.ply ground_mesh.ply - extracts tree locations to file, using a supplied ground mesh" << std::endl;
  std::cout << "                         --tree_roundness 2 - 1: willow, 0.5: birch, 0.2: pine (length per crown radius)." << std::endl;
  std::cout << std::endl;
  std::cout << "rayextract terrain cloud.ply                - extract terrain undersurface to mesh. Slow, so consider decimating first." << std::endl;
  std::cout << std::endl;
  std::cout << "                                 --verbose  - extra debug output." << std::endl;

  exit(error);
}


// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
// #define TEST_CLUSTER
#if defined TEST_CLUSTER
  ray::DebugDraw::init(argc, argv, "rayextract");
  for (int j = 0; j<100; j++)
  {
    std::vector<Eigen::Vector3d> points;
    int num = 1 + std::rand()%4;
    for (int i = 0; i<num; i++)
    {
      Eigen::Vector3d centre((double)(std::rand()%1000)/1000.0, (double)(std::rand()%1000)/1000.0, (double)(std::rand()%1000)/1000.0);
      int count = 1 + std::rand()%30;
      for (int i = 0; i<count; i++)
      {
        Eigen::Vector3d offset((double)(std::rand()%1000)/3000.0, (double)(std::rand()%1000)/3000.0, (double)(std::rand()%1000)/3000.0);
        points.push_back(centre + offset);
      }
    }
    std::cout << "points size: " << points.size() << std::endl;
    ray::generateClusters(points, 10.0, 0.35, true, true);
  }
  return 0;
#endif

  ray::FileArgument cloud_file, mesh_file, trunks_file;
  ray::TextArgument forest("forest"), trees("trees"), trunks("trunks"), terrain("terrain");
  ray::DoubleArgument tree_roundness(0.01, 3.0);
  ray::OptionalKeyValueArgument roundness_option("tree_roundness", 't', &tree_roundness);
  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_trunks = ray::parseCommandLine(argc, argv, {&trunks, &cloud_file}, {&verbose});
  bool extract_trees = ray::parseCommandLine(argc, argv, {&trees, &cloud_file, &trunks_file}, {&verbose});
  bool extract_forest = ray::parseCommandLine(argc, argv, {&forest, &cloud_file, &mesh_file}, {&roundness_option, &verbose});
  bool extract_terrain = ray::parseCommandLine(argc, argv, {&terrain, &cloud_file}, {&verbose});
  if (!extract_trunks && !extract_forest && !extract_terrain && !extract_trees)
    usage();  
  if (verbose.isSet())
  {
    ray::DebugDraw::init(argc, argv, "rayextract");
  }

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage(true);

  if (extract_trunks)
  {
    const double radius = 0.15; // ~ /2 up to *2. So tree diameters 15 cm up to 60 cm 
    ray::Wood woods(cloud, radius, verbose.isSet());
    woods.save(cloud_file.nameStub() + "_trunks.txt");
  }
  else if (extract_trees)
  {
    std::vector<std::pair<Eigen::Vector3d, double> > trunks = ray::Wood::load(trunks_file.name());
    if (trunks.empty())
    {
      std::cerr << "no trunks found in file: " << trunks_file.name() << std::endl;
      usage(true);
    }
    ray::Trees trees(cloud, trunks, verbose.isSet());
    trees.save(cloud_file.nameStub() + "_trees.txt");
  }
  else if (extract_forest)
  {
    ray::Forest forest;
    forest.tree_roundness = roundness_option.isSet() ? tree_roundness.value() : 0.5;
    std::cout << "tree roundness: " << forest.tree_roundness << std::endl;
    forest.verbose = verbose.isSet();

//#define TEST
#if defined TEST
    const int res = 256;
    const double max_tree_height = (double)res / 8.0;
    
    Eigen::ArrayXXd heightfield = Eigen::ArrayXXd::Constant(res, res, -1e10);
    // now lets give it a base hilly floor
    for (int i = 0; i<res; i++)
    {
      for (int j = 0; j<res; j++)
      {
        double h = std::sin(0.1 * ((double)i + 2.0*double(j)));
        heightfield(i,j) = h + ray::random(-0.25, 0.25);
      }
    }

    int num = 250; // 500
    const double radius_to_height = 0.3;//4; // actual radius to height is 2 * radius_to_height^2
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
            double height = p[2] - 1.0/(0.15*p[2]) * mag2; // (p[2]/2.0)*mag2/(radius*radius);
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
    forest.extract(heightfield, mesh, 1.0);
#else
    ray::Mesh mesh;
    ray::readPlyMesh(mesh_file.name(), mesh);
    double voxel_width = 0.25; // 4.0 * cloud.estimatePointSpacing();
    forest.extract(cloud, mesh, voxel_width);
    forest.save(cloud_file.nameStub() + "_trunks.txt");
#endif
  }
  else if (extract_terrain)
  {
    ray::Terrain terrain;
    const double gradient = 1.0; // a half-way divide between ground and wall
    terrain.extract(cloud, cloud_file.nameStub(), gradient, verbose.isSet());
  }
  else
    usage(true);
}

