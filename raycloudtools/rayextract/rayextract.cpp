// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/extraction/raytrunks.h"
#include "raylib/extraction/raybranches.h"
#include "raylib/extraction/rayterrain.h"
#include "raylib/extraction/rayforest.h"
#include "raylib/extraction/raytrees.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayparse.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/extraction/rayclusters.h"
#include "raylib/rayforestgen.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(bool error=false)
{
  std::cout << "Extract natural features into a text file structure" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayextract terrain cloud.ply                - extract terrain undersurface to mesh. Slow, so consider decimating first." << std::endl;
  std::cout << "rayextract trunks cloud.ply                 - extract tree trunk base locations and radii to text file" << std::endl;
  std::cout << "rayextract forest cloud.ply                 - extracts tree locations, radii and heights to file" << std::endl;
  std::cout << "                            --ground ground_mesh.ply - ground mesh file (otherwise assume flat)" << std::endl; 
  std::cout << "                            --trunks cloud_trunks.txt - known tree trunks file" << std::endl;
  std::cout << "                            --width 0.25    - grid cell width" << std::endl;
  std::cout << "                            --smooth 15     - canopy smooth iterations, higher for rough canopies" << std::endl;
  std::cout << "                            --drop_ratio 0.1- here a drop of 10% in canopy height is classed as separate trees" << std::endl;
  std::cout << "rayextract forest_ag cloud.ply              - agglomeration clustering version of above forest extraction method" << std::endl;
  std::cout << "                            --ground ground_mesh.ply - ground mesh file (otherwise assume flat)" << std::endl; 
  std::cout << "                            --trunks cloud_trunks.txt - known tree trunks file" << std::endl;
  std::cout << "                            --width 0.25    - grid cell width" << std::endl;
  std::cout << "                            --min_gradient 0.15 - smallest distance per height to separate clusters" << std::endl;
  std::cout << "                            --max_gradient 1.0  - (-x) largest distance per height to separate clusters" << std::endl;
  std::cout << "rayextract trees cloud.ply cloud_trunks.txt - estimate trees using trunks as seeds, and save to text file" << std::endl;
//  std::cout << "rayextract branches cloud.ply               - estimate tree branches and save to text file" << std::endl;
  std::cout << "                                 --verbose  - extra debug output." << std::endl;

  exit(error);
}


// extracts natural features from scene
int main(int argc, char *argv[])
{ 
  ray::FileArgument cloud_file, mesh_file, trunks_file;
  ray::TextArgument forest("forest"), trees("trees"), trunks("trunks"), branches("branches"), terrain("terrain"), forest_agglomerated("forest_ag");
  ray::OptionalKeyValueArgument groundmesh_option("ground", 'g', &mesh_file);
  ray::OptionalKeyValueArgument trunks_option("trunks", 't', &trunks_file);
  ray::DoubleArgument width(0.01, 10.0), drop(0.001, 1.0), max_gradient(0.01, 5.0), min_gradient(0.01, 5.0);
  ray::IntArgument smooth(0, 50);
  ray::OptionalKeyValueArgument width_option("width", 'w', &width), smooth_option("smooth", 's', &smooth), drop_option("drop_ratio", 'd', &drop);
  ray::OptionalKeyValueArgument max_gradient_option("max_gradient", 'x', &max_gradient), min_gradient_option("min_gradient", 'm', &min_gradient);


  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_terrain = ray::parseCommandLine(argc, argv, {&terrain, &cloud_file}, {&verbose});
  bool extract_trunks = ray::parseCommandLine(argc, argv, {&trunks, &cloud_file}, {&verbose});
  bool extract_forest = ray::parseCommandLine(argc, argv, {&forest, &cloud_file}, {&groundmesh_option, &trunks_option, &width_option, &smooth_option, &drop_option, &verbose});
  bool extract_forest_agglomerate = ray::parseCommandLine(argc, argv, {&forest_agglomerated, &cloud_file}, {&groundmesh_option, &trunks_option, &width_option, &min_gradient_option, &max_gradient_option, &verbose});
  bool extract_trees = ray::parseCommandLine(argc, argv, {&trees, &cloud_file, &trunks_file}, {&verbose});
  bool extract_branches = ray::parseCommandLine(argc, argv, {&branches, &cloud_file}, {&verbose});
  if (!extract_trunks && !extract_branches && !extract_forest && !extract_forest_agglomerate && !extract_terrain && !extract_trees)
    usage();  
  if (verbose.isSet() && (extract_trunks || extract_trees || extract_branches))
  {
    ray::DebugDraw::init(argc, argv, "rayextract");
  }

  if (extract_trunks)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

    const double radius = 0.1; // ~ /2 up to *2. So tree diameters 10 cm up to 40 cm 
    ray::Bush woods(cloud, radius, verbose.isSet(), true);
    woods.save(cloud_file.nameStub() + "_trunks.txt");
  }
/*  else if (extract_branches)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

    const double radius = 0.1; // ~ /2 up to *2. So tree diameters 15 cm up to 60 cm 
    ray::Bush woods(cloud, radius, verbose.isSet(), false);
    woods.save(cloud_file.nameStub() + "_branches.txt");
  }  */
  else if (extract_trees)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

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
    double cell_width = width_option.isSet() ? width.value() : 0.25;
    forest.verbose = verbose.isSet();
    if (smooth_option.isSet())
    {
      forest.smooth_iterations_ = smooth.value();
    }
    if (drop_option.isSet())
    {
      forest.drop_ratio_ = drop.value();
    }
    ray::Mesh mesh;
    if (groundmesh_option.isSet())
    {
      if (!ray::readPlyMesh(mesh_file.name(), mesh))
      {
        std::cerr << "cannot read ground mesh file: " << mesh_file.name() << std::endl;
         usage(true);
      }
    }
    std::vector<std::pair<Eigen::Vector3d, double> > trunks;
    if (trunks_option.isSet())
    {
      trunks = ray::Wood::load(trunks_file.name());
      if (trunks.empty())
      {
        std::cerr << "no trunks found in file: " << trunks_file.name() << std::endl;
        usage(true);
      }
    }
    ray::ForestStructure results = forest.extract(cloud_file.nameStub(), mesh, trunks, cell_width);
    results.save(cloud_file.nameStub() + "_forest.txt");
  }
  else if (extract_forest_agglomerate)
  {
    ray::Forest forest;
    double cell_width = width_option.isSet() ? width.value() : 0.25;
    forest.verbose = verbose.isSet();
    forest.agglomerate_ = true;
    if (max_gradient_option.isSet())
    {
      forest.max_diameter_per_height_ = max_gradient.value();
    }
    if (min_gradient_option.isSet())
    {
      forest.min_diameter_per_height_ = min_gradient.value();     
    }
    ray::Mesh mesh;
    if (groundmesh_option.isSet())
    {
      if (!ray::readPlyMesh(mesh_file.name(), mesh))
      {
        std::cerr << "cannot read ground mesh file: " << mesh_file.name() << std::endl;
         usage(true);
      }
    }
    std::vector<std::pair<Eigen::Vector3d, double> > trunks;
    if (trunks_option.isSet())
    {
      trunks = ray::Wood::load(trunks_file.name());
      if (trunks.empty())
      {
        std::cerr << "no trunks found in file: " << trunks_file.name() << std::endl;
        usage(true);
      }
    }
    ray::ForestStructure results = forest.extract(cloud_file.nameStub(), mesh, trunks, cell_width);
    results.save(cloud_file.nameStub() + "_forest_ag.txt");
  }
  else if (extract_terrain)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

    ray::Terrain terrain;
    const double gradient = 1.0; // a half-way divide between ground and wall
    terrain.extract(cloud, cloud_file.nameStub(), gradient, verbose.isSet());
  }
  else
    usage(true);
}

