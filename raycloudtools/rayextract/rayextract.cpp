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
  std::cout << "                            --gradient 1    - maximum gradient counted as terrain" << std::endl;
  std::cout << "rayextract trunks cloud.ply                 - extract tree trunk base locations and radii to text file" << std::endl;
  std::cout << "                            --exclude_rays  - does not use rays to exclude candidates with rays passing through" << std::endl;
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
  std::cout << "rayextract trees cloud.ply ground_mesh.ply  - estimate trees, and save to text file" << std::endl;
  std::cout << "                            --max_diameter 0.9   - (-m) maximum trunk diameter in segmenting trees" << std::endl;
  std::cout << "                            --min_diameter 0.02  - (-n) minimum branch diameter" << std::endl;
  std::cout << "                            --distance_limit 1   - (-d) maximum distance between neighbour points in a tree" << std::endl;
  std::cout << "                            --height_min 2       - (-h) minimum height counted as a tree" << std::endl;
  std::cout << "                            --min_length_per_radius 140- (-l) the tapering rate of branches" << std::endl;
//  std::cout << "                            --radius_exponent 0.67 - (-e) exponent of radius in estimating length" << std::endl;
//  std::cout << "                            --cylinder_length_to_width 4- (-c) how slender the cylinders are" << std::endl;
//  std::cout << "                            --gap_ratio 2.5      - (-g) will split for lateral gaps at this multiple of radius" << std::endl;
//  std::cout << "                            --span_ratio 4.5     - (-s) will split when branch width spans this multiple of radius" << std::endl;
  std::cout << "                            --gravity_factor 0.3 - (-f) larger values preference vertical trees" << std::endl;
  std::cout << "                            --branch_segmentation- (-b) _segmented.ply is per branch segment" << std::endl;
  std::cout << "                            --grid_width         - (-w) crops results assuming cloud has been gridded with given width" << std::endl;

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
  ray::DoubleArgument gradient(0.001, 1000.0);
  ray::OptionalKeyValueArgument gradient_option("gradient", 'g', &gradient);
  ray::OptionalFlagArgument exclude_rays("exclude_rays", 'e'), segment_branches("branch_segmentation", 'b');
  ray::DoubleArgument width(0.01, 10.0), drop(0.001, 1.0), max_gradient(0.01, 5.0), min_gradient(0.01, 5.0);

  ray::DoubleArgument max_diameter(0.01, 100.0), distance_limit(0.01, 10.0), height_min(0.01, 1000.0), min_diameter(0.01, 100.0);
  ray::DoubleArgument length_to_radius(0.01, 10000.0), cylinder_length_to_width(0.1, 20.0), gap_ratio(0.01, 10.0), span_ratio(0.01, 10.0);
  ray::DoubleArgument gravity_factor(0.0, 100.0), radius_exponent(0.0, 100.0), grid_width(1.0, 100000.0), grid_overlap(0.0, 0.9);
  ray::OptionalKeyValueArgument max_diameter_option("max_diameter", 'm', &max_diameter);
  ray::OptionalKeyValueArgument min_diameter_option("min_diameter", 'n', &min_diameter);
  ray::OptionalKeyValueArgument distance_limit_option("distance_limit", 'd', &distance_limit);
  ray::OptionalKeyValueArgument height_min_option("height_min", 'h', &height_min);
  ray::OptionalKeyValueArgument length_to_radius_option("min_length_per_radius", 'l', &length_to_radius);
  ray::OptionalKeyValueArgument radius_exponent_option("radius_exponent", 'e', &radius_exponent);
  ray::OptionalKeyValueArgument cylinder_length_to_width_option("cylinder_length_to_width", 'c', &cylinder_length_to_width);
  ray::OptionalKeyValueArgument gap_ratio_option("gap_ratio", 'g', &gap_ratio);
  ray::OptionalKeyValueArgument span_ratio_option("span_ratio", 's', &span_ratio);
  ray::OptionalKeyValueArgument gravity_factor_option("gravity_factor", 'f', &gravity_factor);
  ray::OptionalKeyValueArgument grid_width_option("grid_width", 'w', &grid_width);

  ray::IntArgument smooth(0, 50);
  ray::OptionalKeyValueArgument width_option("width", 'w', &width), smooth_option("smooth", 's', &smooth), drop_option("drop_ratio", 'd', &drop);
  ray::OptionalKeyValueArgument max_gradient_option("max_gradient", 'x', &max_gradient), min_gradient_option("min_gradient", 'm', &min_gradient);


  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_terrain = ray::parseCommandLine(argc, argv, {&terrain, &cloud_file}, {&gradient_option, &verbose});
  bool extract_trunks = ray::parseCommandLine(argc, argv, {&trunks, &cloud_file}, {&exclude_rays, &verbose});
  bool extract_forest = ray::parseCommandLine(argc, argv, {&forest, &cloud_file}, {&groundmesh_option, &trunks_option, &width_option, &smooth_option, &drop_option, &verbose});
  bool extract_forest_agglomerate = ray::parseCommandLine(argc, argv, {&forest_agglomerated, &cloud_file}, {&groundmesh_option, &trunks_option, &width_option, &min_gradient_option, &max_gradient_option, &verbose});
  bool extract_trees = ray::parseCommandLine(argc, argv, {&trees, &cloud_file, &mesh_file}, {&max_diameter_option, &distance_limit_option, &height_min_option, &min_diameter_option, &length_to_radius_option, &cylinder_length_to_width_option, &gap_ratio_option, &span_ratio_option, &gravity_factor_option, &radius_exponent_option, &segment_branches, &grid_width_option, &verbose});
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
    ray::Bush woods(cloud, radius, verbose.isSet(), true, exclude_rays.isSet());
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

    ray::Mesh mesh;
    if (!ray::readPlyMesh(mesh_file.name(), mesh))
      usage(true);
    ray::TreesParams params;
    if (max_diameter_option.isSet())
      params.max_diameter = max_diameter.value();
    if (distance_limit_option.isSet())
      params.distance_limit = distance_limit.value();
    if (height_min_option.isSet())
      params.height_min = height_min.value();
    if (min_diameter_option.isSet())
      params.min_diameter = min_diameter.value();
    if (length_to_radius_option.isSet())
      params.length_to_radius = length_to_radius.value();
    if (radius_exponent_option.isSet())
      params.radius_exponent = radius_exponent.value();
    if (cylinder_length_to_width_option.isSet())
      params.cylinder_length_to_width = cylinder_length_to_width.value();
    if (gap_ratio_option.isSet())
      params.gap_ratio = gap_ratio.value();
    if (span_ratio_option.isSet())
      params.span_ratio = span_ratio.value();    
    if (gravity_factor_option.isSet())
      params.gravity_factor = gravity_factor.value();
    if (grid_width_option.isSet())
      params.grid_width = grid_width.value();
    params.segment_branches = segment_branches.isSet();
  
    ray::Trees trees(cloud, mesh, params, verbose.isSet());

    trees.save(cloud_file.nameStub() + "_trees.txt");
    cloud.save(cloud_file.nameStub() + "_segmented.ply");
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
        usage(true);
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
    const double grad = gradient_option.isSet() ? gradient.value() : 1.0; 
    terrain.extract(cloud, cloud_file.nameStub(), grad, verbose.isSet());
  }
  else
    usage(true);
}

