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
  std::cout << "rayextract trees cloud.ply cloud_trunks.txt - estimate trees using trunks as seeds, and save to text file" << std::endl;
  std::cout << "rayextract branches cloud.ply               - estimate tree branches and save to text file" << std::endl;
  std::cout << "                                 --verbose  - extra debug output." << std::endl;

  exit(error);
}


// extracts natural features from scene
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file, mesh_file, trunks_file;
  ray::TextArgument forest("forest"), trees("trees"), trunks("trunks"), branches("branches"), terrain("terrain");
  ray::OptionalKeyValueArgument groundmesh_option("ground", 'g', &mesh_file);
  ray::OptionalKeyValueArgument trunks_option("trunks", 't', &trunks_file);
  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_terrain = ray::parseCommandLine(argc, argv, {&terrain, &cloud_file}, {&verbose});
  bool extract_trunks = ray::parseCommandLine(argc, argv, {&trunks, &cloud_file}, {&verbose});
  bool extract_forest = ray::parseCommandLine(argc, argv, {&forest, &cloud_file}, {&groundmesh_option, &trunks_option, &verbose});
  bool extract_trees = ray::parseCommandLine(argc, argv, {&trees, &cloud_file, &trunks_file}, {&verbose});
  bool extract_branches = ray::parseCommandLine(argc, argv, {&branches, &cloud_file}, {&verbose});
  if (!extract_trunks && !extract_branches && !extract_forest && !extract_terrain && !extract_trees)
    usage();  
  if (verbose.isSet())
  {
    ray::DebugDraw::init(argc, argv, "rayextract");
  }


  if (extract_trunks)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

 /*   const double radius = 0.15; // ~ /2 up to *2. So tree diameters 15 cm up to 60 cm 
    ray::Wood woods(cloud, radius, verbose.isSet());*/

    const double radius = 0.1; // ~ /2 up to *2. So tree diameters 10 cm up to 40 cm 
    ray::Bush woods(cloud, radius, verbose.isSet(), true);
    woods.save(cloud_file.nameStub() + "_trunks.txt");
  }
  else if (extract_branches)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
      usage(true);

    const double radius = 0.1; // ~ /2 up to *2. So tree diameters 15 cm up to 60 cm 
    ray::Bush woods(cloud, radius, verbose.isSet(), false);
    woods.save(cloud_file.nameStub() + "_branches.txt");
  }  
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
    forest.verbose = verbose.isSet();

    ray::Mesh mesh;
    if (groundmesh_option.isSet())
    {
      if (!ray::readPlyMesh(mesh_file.name(), mesh))
      {
        std::cerr << "cannot read ground mesh file: " << mesh_file.name() << std::endl;
         usage(true);
      }
    }
    std::vector<std::pair<Eigen::Vector3d, double> > trunks = ray::Wood::load(trunks_file.name());
    if (trunks.empty())
    {
      std::cerr << "no trunks found in file: " << trunks_file.name() << std::endl;
      usage(true);
    }

    std::vector<ray::TreeSummary> results = forest.extract(cloud_file.nameStub(), mesh, trunks);
    ray::TreeSummary::save(cloud_file.nameStub() + "_forest.txt", results);
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

