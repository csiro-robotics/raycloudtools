// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/extraction/rayclusters.h"
#include "raylib/extraction/rayforest.h"
#include "raylib/extraction/rayterrain.h"
#include "raylib/extraction/raytrees.h"
#include "raylib/extraction/raytrunks.h"
#include "raylib/extraction/rayleaves.h"
#include "raylib/raycloud.h"
#include "raylib/rayforestgen.h"
#include "raylib/rayforeststructure.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

static std::string extract_type;

void usage(int exit_code = 1)
{
  const bool none = extract_type != "terrain" && extract_type != "trunks" && extract_type != "forest" && extract_type != "trees" && extract_type != "leaves";
  // clang-format off
  std::cout << "Extract natural features into a text file or mesh file" << std::endl;
  std::cout << "usage:" << std::endl;
  if (extract_type == "terrain" || none)
  {
    std::cout << "rayextract terrain cloud.ply                - extract terrain undersurface to mesh. Slow, so consider decimating first." << std::endl;
    std::cout << "                            --gradient 1    - maximum gradient counted as terrain" << std::endl;
  }
  if (extract_type == "trunks" || none)
  {
    std::cout << "rayextract trunks cloud.ply                 - extract tree trunk base locations and radii to text file" << std::endl;
    std::cout << "                            --exclude_rays  - does not use rays to exclude candidates with rays passing through" << std::endl;
  }
  if (extract_type == "forest" || none)
  {
    std::cout << "rayextract forest cloud.ply                 - extracts tree locations, radii and heights to file" << std::endl;
    std::cout << "                            --ground ground_mesh.ply - ground mesh file (otherwise assume flat)" << std::endl; 
    std::cout << "                            --trunks cloud_trunks.txt - known tree trunks file" << std::endl;
    std::cout << "                            --width 0.25    - grid cell width" << std::endl;
    std::cout << "                            --smooth 15     - canopy smooth iterations, higher for rough canopies" << std::endl;
    std::cout << "                            --drop_ratio 0.1- here a drop of 10% in canopy height is classed as separate trees" << std::endl;
  }
  if (extract_type == "trees" || none)
  {
    std::cout << "rayextract trees cloud.ply ground_mesh.ply  - estimate trees, and save to text file, mesh file, and segmented (coloured per-tree) cloud. Use TreeTools (github) to manipulate text file." << std::endl;
    std::cout << "                            --max_diameter 0.9   - (-m) maximum trunk diameter in segmenting trees" << std::endl;
    std::cout << "                            --crop_length 1.0    - (-p) crops small branches to this distance from end" << std::endl;
    std::cout << "                            --distance_limit 1   - (-d) maximum distance between neighbour points in a tree" << std::endl;
    std::cout << "                            --height_min 2       - (-h) minimum height counted as a tree" << std::endl;
    std::cout << "                            --girth_height_ratio 0.12 - (-i) the amount up tree's height to estimate trunk girth" << std::endl;
    std::cout << "                            --global_taper 0.024 - (-a) force a taper value (diameter per length) for trees under global_taper_factor of max tree height. Use 0 to estimate global taper from the data" << std::endl;
    std::cout << "                            --global_taper_factor 0.3- (-o) 1 estimates same taper for whole scan, 0 is per-tree tapering. Like a soft cutoff at this amount of max tree height" << std::endl;
    std::cout << "                            --gravity_factor 0.3 - (-f) larger values preference vertical trees" << std::endl;
    std::cout << "                            --branch_segmentation- (-b) _segmented.ply is per branch segment" << std::endl;
    std::cout << "                            --grid_width 10      - (-w) crops results assuming cloud has been gridded with given width" << std::endl;
    std::cout << "                            --use_rays           - (-u) use rays to reduce trunk radius overestimation in noisy cloud data" << std::endl;
    std::cout << "                            (for internal constants -c -g -s see source file rayextract)" << std::endl;
  // These are the internal parameters that I don't expose as they are 'advanced' only, you shouldn't need to adjust them
  //  std::cout << "                            --cylinder_length_to_width 4- (-c) how slender the cylinders are" << std::endl;
  //  std::cout << "                            --gap_ratio 0.016    - (-g) will split for lateral gaps at this multiple of branch length" << std::endl;
  //  std::cout << "                            --span_ratio 4.5     - (-s) will split when branch width spans this multiple of radius" << std::endl;
  }
  if (extract_type == "leaves" || none)
  {
    std::cout << "rayextract leaves cloud.ply trees.txt            - reconstruct the leaf locations coming from the specified tree structures, and save to text file" << std::endl;
    std::cout << "                            --leaf mesh.ply      - mesh for each leaf. Should have its centre at 0,0,0, be along the y axis, with first vertex at the stalk connection" << std::endl;
    std::cout << "                            --leaf image.png     - make leaf from the image, assuming it uses its alpha channel (loads in cloudompare but not meshlab)" << std::endl;
    std::cout << "                            --leaf_area 0.002    - area for each leaf." << std::endl;
    std::cout << "                            --leaf_droop 0.1     - drop per square horizontal distance." << std::endl;
    std::cout << "                            --stalks             - include stalks to closest branch." << std::endl;
    std::cout << "                                 --verbose  - extra debug output." << std::endl;
  }
  // clang-format on
  exit(exit_code);
}


/// extracts natural features from a scene
int rayExtract(int argc, char *argv[])
{
  if (argc > 1)
  {
    extract_type = std::string(argv[1]);
  }
  ray::FileArgument cloud_file, mesh_file, trunks_file, trees_file, leaf_file;
  ray::TextArgument forest("forest"), trees("trees"), trunks("trunks"), terrain("terrain"), leaves("leaves");
  ray::OptionalKeyValueArgument groundmesh_option("ground", 'g', &mesh_file);
  ray::OptionalKeyValueArgument trunks_option("trunks", 't', &trunks_file);
  ray::DoubleArgument gradient(0.001, 1000.0, 1.0), global_taper(0.0, 1.0), global_taper_factor(0.0, 1.0);
  ray::OptionalKeyValueArgument gradient_option("gradient", 'g', &gradient);
  ray::OptionalFlagArgument exclude_rays("exclude_rays", 'e'), segment_branches("branch_segmentation", 'b'), stalks("stalks", 's'), use_rays("use_rays", 'u');
  ray::DoubleArgument width(0.01, 10.0, 0.25), drop(0.001, 1.0), max_gradient(0.01, 5.0), min_gradient(0.01, 5.0);

  ray::DoubleArgument max_diameter(0.01, 100.0), distance_limit(0.01, 10.0), height_min(0.01, 1000.0),
    min_diameter(0.01, 100.0), leaf_area(0.00001, 1.0, 0.002), leaf_droop(0.0, 10.0, 0.1), crop_length(0.01, 100.0);;
  ray::DoubleArgument girth_height_ratio(0.001, 0.5), length_to_radius(0.01, 10000.0), cylinder_length_to_width(0.1, 20.0), gap_ratio(0.01, 10.0),
    span_ratio(0.01, 10.0);
  ray::DoubleArgument gravity_factor(0.0, 100.0), grid_width(1.0, 100000.0),
    grid_overlap(0.0, 0.9);
  ray::OptionalKeyValueArgument max_diameter_option("max_diameter", 'm', &max_diameter);
  ray::OptionalKeyValueArgument crop_length_option("crop_length", 'n', &crop_length);
  ray::OptionalKeyValueArgument distance_limit_option("distance_limit", 'd', &distance_limit);
  ray::OptionalKeyValueArgument height_min_option("height_min", 'h', &height_min);
  ray::OptionalKeyValueArgument girth_height_ratio_option("girth_height_ratio", 'i', &girth_height_ratio);
  ray::OptionalKeyValueArgument cylinder_length_to_width_option("cylinder_length_to_width", 'c',
                                                                &cylinder_length_to_width);
  ray::OptionalKeyValueArgument gap_ratio_option("gap_ratio", 'g', &gap_ratio);
  ray::OptionalKeyValueArgument span_ratio_option("span_ratio", 's', &span_ratio);
  ray::OptionalKeyValueArgument gravity_factor_option("gravity_factor", 'f', &gravity_factor);
  ray::OptionalKeyValueArgument grid_width_option("grid_width", 'w', &grid_width);
  ray::OptionalKeyValueArgument global_taper_option("global_taper", 'a', &global_taper);
  ray::OptionalKeyValueArgument global_taper_factor_option("global_taper_factor", 'o', &global_taper_factor);
  ray::OptionalKeyValueArgument leaf_option("leaf", 'l', &leaf_file);
  ray::OptionalKeyValueArgument leaf_area_option("leaf_area", 'a', &leaf_area);
  ray::OptionalKeyValueArgument leaf_droop_option("leaf_droop", 'd', &leaf_droop);

  ray::IntArgument smooth(0, 50);
  ray::OptionalKeyValueArgument width_option("width", 'w', &width), smooth_option("smooth", 's', &smooth),
    drop_option("drop_ratio", 'd', &drop);

  ray::OptionalFlagArgument verbose("verbose", 'v');

  bool extract_terrain = ray::parseCommandLine(argc, argv, { &terrain, &cloud_file }, { &gradient_option, &verbose });
  bool extract_trunks = ray::parseCommandLine(argc, argv, { &trunks, &cloud_file }, { &exclude_rays, &verbose });
  bool extract_forest = ray::parseCommandLine(
    argc, argv, { &forest, &cloud_file },
    { &groundmesh_option, &trunks_option, &width_option, &smooth_option, &drop_option, &verbose });
  bool extract_trees = ray::parseCommandLine(
    argc, argv, { &trees, &cloud_file, &mesh_file },
    { &max_diameter_option, &distance_limit_option, &height_min_option, &crop_length_option, &girth_height_ratio_option,
      &cylinder_length_to_width_option, &gap_ratio_option, &span_ratio_option, &gravity_factor_option,
      &segment_branches, &grid_width_option, &global_taper_option, &global_taper_factor_option, &use_rays, &verbose });
  bool extract_leaves = ray::parseCommandLine(argc, argv, { &leaves, &cloud_file, &trees_file }, { &leaf_option, &leaf_area_option, &leaf_droop_option, &stalks });


  if (!extract_trunks && !extract_forest && !extract_terrain && !extract_trees && !extract_leaves)
  {
    usage();
  }

  // finds cylindrical trunks in the data and saves them to an _trunks.txt file
  if (extract_trunks)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
    {
      usage(true);
    }
    Eigen::Vector3d offset = cloud.removeStartPos();

    const double radius = 0.1;  // ~ /2 up to *2. So tree diameters 10 cm up to 40 cm
    ray::Trunks trunks(cloud, offset, radius, verbose.isSet(), exclude_rays.isSet());
    trunks.save(cloud_file.nameStub() + "_trunks.txt", offset);
  }
  // finds full tree structures (piecewise cylindrical representation) and saves to file
  else if (extract_trees)
  {
    ray::Cloud cloud;
    const int min_num_rays = 40;
    if (!cloud.load(cloud_file.name(), true, min_num_rays))
    {
      usage(true);
    }
    Eigen::Vector3d offset = cloud.removeStartPos();

    ray::Mesh mesh;
    if (!ray::readPlyMesh(mesh_file.name(), mesh))
    {
      usage(true);
    }
    mesh.translate(-offset);

    ray::TreesParams params;
    if (max_diameter_option.isSet())
    {
      params.max_diameter = max_diameter.value();
    }
    if (distance_limit_option.isSet())
    {
      params.distance_limit = distance_limit.value();
    }
    if (height_min_option.isSet())
    {
      params.height_min = height_min.value();
    }
    if (crop_length_option.isSet())
    {
      params.crop_length = crop_length.value();
    }
    if (girth_height_ratio_option.isSet())
    {
      params.girth_height_ratio = girth_height_ratio.value();
    }
    if (cylinder_length_to_width_option.isSet())
    {
      params.cylinder_length_to_width = cylinder_length_to_width.value();
    }
    if (gap_ratio_option.isSet())
    {
      params.gap_ratio = gap_ratio.value();
    }
    if (span_ratio_option.isSet())
    {
      params.span_ratio = span_ratio.value();
    }
    if (gravity_factor_option.isSet())
    {
      params.gravity_factor = gravity_factor.value();
    }
    if (grid_width_option.isSet())
    {
      params.grid_width = grid_width.value();
    }
    if (global_taper_option.isSet())
    {
      params.global_taper = 0.5 * global_taper.value();
    }
    if (global_taper_factor_option.isSet())
    {
      params.global_taper_factor = global_taper_factor.value();
    }   
    params.use_rays = use_rays.isSet(); 
    params.segment_branches = segment_branches.isSet();

    ray::Trees trees(cloud, offset, mesh, params, verbose.isSet());

    // output the picewise cylindrical description of the trees
    trees.save(cloud_file.nameStub() + "_trees.txt", offset, verbose.isSet());
    // we also save a segmented (one colour per tree) file, as this is a useful output
    cloud.translate(offset);
    cloud.save(cloud_file.nameStub() + "_segmented.ply");
    // let's also save the trees out as a mesh
    // it is a bit inefficient to load from file just to convert it into the forest structure, but
    // it works OK for now. Better would be for ray::Trees so store the result as a ray::ForestStructure
    ray::ForestStructure forest;
    if (!forest.load(cloud_file.nameStub() + "_trees.txt"))
    {
      usage();
    }
    ray::Mesh tree_mesh;
    forest.generateSmoothMesh(tree_mesh, -1, 1, 1, 1);
    ray::writePlyMesh(cloud_file.nameStub() + "_trees_mesh.ply", tree_mesh, true);    
  }
  // extract the tree locations from a larger, aerial view of a forest
  else if (extract_forest)
  {
    ray::Forest forest;
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
        usage(true);
      }
    }
    std::vector<std::pair<Eigen::Vector3d, double>> trunks;
    // the results from extracting trunks can optionally be passed in, as a guide
    if (trunks_option.isSet())
    {
      ray::ForestStructure forest;
      if (!forest.load(trunks_file.name()))
      {
        usage(true);
      }
      for (auto &tree : forest.trees)
      {
        trunks.push_back(std::pair<Eigen::Vector3d, double>(tree.segments()[0].tip, tree.segments()[0].radius));
      }
    }
    ray::ForestStructure results = forest.extract(cloud_file.nameStub(), mesh, trunks, width.value());
    // save the results, which is a location, radius and height per tree
    results.save(cloud_file.nameStub() + "_forest.txt");
  }
  // extract the terrain to a .ply mesh file
  // this uses a sand model (no terrain is sloped more than 'gradient') which is a
  // highest lower bound
  else if (extract_terrain)
  {
    ray::Cloud cloud;
    if (!cloud.load(cloud_file.name()))
    {
      usage(true);
    }
    Eigen::Vector3d offset = cloud.removeStartPos();

    ray::Terrain terrain;
    terrain.extract(cloud, offset, cloud_file.nameStub(), gradient.value(), verbose.isSet());
  }
  else if (extract_leaves)
  {
    ray::generateLeaves(cloud_file.nameStub(), trees_file.name(), leaf_file.name(), 
      leaf_area.value(), leaf_droop.value(), stalks.isSet());
  }
  else
  {
    usage(true);
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayExtract, argc, argv);
}