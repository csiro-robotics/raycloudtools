// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "raylib/raycloud.h"
#include "raylib/raycuboid.h"
#include "raylib/raylibconfig.h"
#include "raylib/rayparse.h"
#include "raylib/rayrenderer.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Render a ray cloud as an image, from a specified viewpoint" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrender raycloudfile.ply top ends        - render from the top (plan view) the end points" << std::endl;
  std::cout << "                           left            - facing negative x axis" << std::endl;
  std::cout << "                           right           - facing positive x axis" << std::endl;
  std::cout << "                           front           - facing negative y axis" << std::endl;
  std::cout << "                           back            - facing positive y axis" << std::endl;
  std::cout << "                               mean        - mean colour on axis" << std::endl;
  std::cout << "                               sum         - sum colours (globally scaled to colour range)" << std::endl;
  std::cout << "                               starts      - render the ray start points" << std::endl;
  std::cout << "                               rays        - render the full set of rays" << std::endl;
  std::cout << "                               height      - render the maximum heights in the view axis" << std::endl;
  std::cout << "                               density     - shade according to estimated density within pixel" << std::endl;
  std::cout << "                               density_rgb - r->g->b colour by estimated density" << std::endl;
  std::cout << "                                           - Note: densities weight by return intensity (alpha) with alpha=100 for non-returns." << std::endl;
  std::cout << "                                           - use rayimport --max_intensity to scale, or raycolour alpha 1 to disable weighting." << std::endl;
  std::cout << "                     --resolution 512      - long axis resolution" << std::endl;
  std::cout << "                     --output name.png     - optional output file name. " << std::endl;
  std::cout << "                                             Supports .png, .tga, .hdr, .jpg, .bmp" << std::endl;
  std::cout << "                     --mark_origin         - place a 255,0,255 pixel at the coordinate origin. " << std::endl;
  std::cout << "                     --output_transform    - generate a yaml file containing the" << std::endl;
  std::cout << "                                             transform from the raycloud to" << std::endl;
  std::cout << "                                             pixels. Only compatible with top" << std::endl;
  std::cout << "                                             view." << std::endl;
  std::cout << "                     --georeference name.proj- projection file name, to output (geo)tif file. " << std::endl;
  std::cout << "                     --pixel_width 0.1     - optional pixel width in m, instead of resolution" << std::endl;
  std::cout << "                     --grid_width 100      - optionally bound to a grid cell width such that one cell centre is 0,0" << std::endl;
  std::cout << "Default output is raycloudfile.png" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayRender(int argc, char *argv[])
{
  ray::KeyChoice viewpoint({ "top", "left", "right", "front", "back" });
  ray::KeyChoice style({ "ends", "mean", "sum", "starts", "rays", "height", "density", "density_rgb" });
  ray::DoubleArgument pixel_width(0.0001, 1000.0), grid_width(0.01, 1000000.0);
  ray::IntArgument resolution(1,20000, 512);
  ray::FileArgument cloud_file, image_file, transform_file, projection_file(false);
  ray::OptionalFlagArgument mark_origin("mark_origin", 'm');
  ray::OptionalKeyValueArgument resolution_option("resolution", 'r', &resolution);
  ray::OptionalKeyValueArgument pixel_width_option("pixel_width", 'p', &pixel_width);
  ray::OptionalKeyValueArgument grid_width_option("grid_width", 'g', &grid_width);
  ray::OptionalKeyValueArgument output_file_option("output", 'o', &image_file);
  ray::OptionalKeyValueArgument projection_file_option("georeference", 'g', &projection_file);
  ray::OptionalKeyValueArgument transform_file_option("output_transform", 't', &transform_file);
  if (!ray::parseCommandLine(
        argc, argv, { &cloud_file, &viewpoint, &style },
        { &resolution_option, &pixel_width_option, &output_file_option, &mark_origin, &transform_file_option, &grid_width_option, &projection_file_option }))
  {
    usage();
  }
  if (!output_file_option.isSet())
  {
    image_file.name() = cloud_file.nameStub() + (projection_file_option.isSet() ? ".tif" : ".png");
  }
  // a projection file describes where the ray cloud is in the world, which allows
  // images to be output in geotiff (geolocalised tiff) format.
  if (projection_file_option.isSet())
  {
#if !RAYLIB_WITH_TIFF
    std::cerr << "Error: georeferencing requires the WITH_TIFF build flag enabled. See README.md." << std::endl;
    usage();
#endif
    if (image_file.nameExt() != "tif")
    {
      std::cerr << "Error: projection files can only be used when outputting a .tif file" << std::endl;
      usage();
    }
    if (viewpoint.selectedKey() != "top")
    {
      std::cerr << "Error: can only geolocate a top-down render" << std::endl;
      usage();
    }
  }

  ray::Cloud::Info info;
  if (!ray::Cloud::getInfo(cloud_file.name(), info))
  {
    usage();
  }
  ray::Cuboid bounds = info.ends_bound;  // exclude the unbounded ray lengths (e.g. up into the sky)
  if (grid_width_option.isSet()) // adjust min_bound to be well-aligned
  {
    Eigen::Vector3d mid = (bounds.min_bound_ + bounds.max_bound_)/2.0;
    Eigen::Vector3d min_bound = bounds.min_bound_;
    Eigen::Vector3d max_bound = bounds.max_bound_;
    min_bound[0] = grid_width.value() * std::round(mid[0] / grid_width.value()) - 0.5*grid_width.value();
    min_bound[1] = grid_width.value() * std::round(mid[1] / grid_width.value()) - 0.5*grid_width.value();
    max_bound[0] = min_bound[0] + grid_width.value();
    max_bound[1] = min_bound[1] + grid_width.value();
    if (min_bound[0] > bounds.min_bound_[0] || min_bound[1] > bounds.min_bound_[1] || 
        max_bound[0] < bounds.max_bound_[0] || max_bound[1] < bounds.max_bound_[1])
    {
      std::cout << "Warning: cloud overlaps grid cell of width: " << grid_width.value() << " image bounds extended" << std::endl;
    }
    bounds.min_bound_ = ray::minVector(bounds.min_bound_, min_bound);
    bounds.max_bound_ = ray::maxVector(bounds.max_bound_, max_bound);
  }  
  double pix_width = pixel_width.value();
  if (!pixel_width_option.isSet())
  {
    Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
    double length = std::max(extent[0], extent[1]);
    pix_width = length / resolution.value();
  }
  if (pix_width <= 0.0)
  {
    usage();
  }

  // quick casting allowed, taking care that the text and enums are in the same order
  const ray::ViewDirection view_dir = static_cast<ray::ViewDirection>(viewpoint.selectedID());
  const ray::RenderStyle render_style = static_cast<ray::RenderStyle>(style.selectedID());

  // an option to output the transformation from image to world frame
  if (transform_file_option.isSet() && (view_dir != ray::ViewDirection::Top))
  {
    std::cout << "--output_transform can only be used when view is top." << std::endl;
    usage();
  }

  if (!ray::renderCloud(cloud_file.name(), bounds, view_dir, render_style, pix_width, image_file.name(),
                        projection_file.name(), mark_origin.isSet(),
                        transform_file_option.isSet() ? &transform_file.name() : nullptr))
  {
    usage();
  }

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayRender, argc, argv);
}