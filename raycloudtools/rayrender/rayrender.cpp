// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/raycloud.h"
#include "raylib/rayrenderer.h"

void usage(int exit_code = 1)
{
  std::cout << "Render a ray cloud as an image, from a specified viewpoint" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrender raycloudfile.ply top ends        - render from the top (plan view) the end points" << std::endl;
  std::cout << "                           left            - facing negative x axis" << std::endl;
  std::cout << "                           right           - facing positive x axis" << std::endl;
  std::cout << "                           front           - facing negative y axis" << std::endl;
  std::cout << "                           back            - facing positive y axis" << std::endl;
  std::cout << "                               mean        - mean colour on axis" << std::endl;
  std::cout << "                               sum         - sum colours (globally scaled to colour range)"
    << std::endl;
  std::cout << "                               starts      - render the ray start points" << std::endl;
  std::cout << "                               rays        - render the full set of rays" << std::endl;
  std::cout << "                               density     - shade according to estimated density within pixel"
    << std::endl;
  std::cout << "                               density_rgb - r->g->b colour by estimated density"
    << std::endl;
  std::cout << "                     --pixel_width 0.1     - optional pixel width in m" << std::endl;
  std::cout << "                     --output name.png     - optional output file name. " << std::endl;
  std::cout << "                                             Supports .png, .tga, .hdr, .jpg, .bmp" << std::endl;
  std::cout << "                     --mark_origin         - place a 255,0,255 pixel at the coordinate origin. " << std::endl;
  std::cout << "                     --output_transform    - generate a yaml file containing the" << std::endl;
  std::cout << "                                             transform from the raycloud to" << std::endl;
  std::cout << "                                             pixels. Only compatible with top" << std::endl;
  std::cout << "                                             view." << std::endl;
  std::cout << "Default output is raycloudfile.png" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::KeyChoice viewpoint({"top", "left", "right", "front", "back"});
  ray::KeyChoice style({"ends", "mean", "sum", "starts", "rays", "density", "density_rgb"});
  ray::DoubleArgument pixel_width(0.0001, 1000.0);
  ray::FileArgument cloud_file, image_file, transform_file;
  ray::OptionalFlagArgument mark_origin("mark_origin", 'm');
  ray::OptionalKeyValueArgument pixel_width_option("pixel_width", 'p', &pixel_width);
  ray::OptionalKeyValueArgument output_file_option("output", 'o', &image_file);
  ray::OptionalKeyValueArgument transform_file_option("output_transform", 't', &transform_file);
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &viewpoint, &style}, {&pixel_width_option, &output_file_option, &mark_origin, &transform_file_option}))
  {
    usage();
  }
  if (!output_file_option.isSet())
  {
    image_file.name() = cloud_file.nameStub() + ".png";
  }

  ray::Cloud::Info info;
  if (!ray::Cloud::getInfo(cloud_file.name(), info))
  {
    usage();
  }
  const ray::Cuboid bounds = info.ends_bound; // exclude the unbounded ray lengths (e.g. up into the sky)
  double pix_width = pixel_width.value();
  if (!pixel_width_option.isSet())
  {
    const double spacing_scale = 2.0; // a reasonable default multiplier on the spacing between points
    pix_width = spacing_scale * ray::Cloud::estimatePointSpacing(cloud_file.name(), bounds, info.num_bounded);
  }
  if (pix_width <= 0.0)
  {
    usage();
  }

  // quick casting allowed, taking care that the text and enums are in the same order
  const ray::ViewDirection view_dir = static_cast<ray::ViewDirection>(viewpoint.selectedID());
  const ray::RenderStyle render_style = static_cast<ray::RenderStyle>(style.selectedID());

  if (transform_file_option.isSet() && (view_dir != ray::ViewDirection::Top))
  {
    std::cout << "--output_transform can only be used when view is top." << std::endl;
    usage();
  }

  if (!ray::renderCloud(cloud_file.name(), bounds, view_dir, render_style, pix_width, image_file.name(),
                        mark_origin.isSet(), transform_file_option.isSet() ? &transform_file.name() : nullptr))
  {
    usage();
  }

  return 0;
}
