// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imagewrite.h"

void usage(int exit_code = 1)
{
  std::cout << "Render a ray cloud as an image, from a specified viewpoint" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayrender raycloudfile top ends outputfile.png - render file from top (plan view), in solid colour" << std::endl;
  std::cout << "                       left                     - facing negative x axis" << std::endl;
  std::cout << "                       right                    - facing positive x axis" << std::endl;
  std::cout << "                       front                    - facing negative y axis" << std::endl;
  std::cout << "                       back                     - facing positive y axis" << std::endl;
  std::cout << "                           mean                 - mean colour on axis" << std::endl;
  std::cout << "                           sum                  - sum colours (globally scaled to colour range)"
    << std::endl;
  std::cout << "                           density              - shade according to estimated density within pixel"
    << std::endl;
  std::cout << "                           starts               - render the ray start points" << std::endl;
  std::cout << "                           rays                 - render the full set of rays" << std::endl;
  std::cout << "                                 outputfile.hdr - format allows a wider scale range"
    << std::endl;
  std::cout << "                              --pixel_width 0.1 - optional pixel width in m" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::KeyChoice viewpoint({"top", "left", "right", "front", "back"});
  ray::KeyChoice style({"ends", "mean", "sum", "density", "starts", "rays"});
  ray::DoubleArgument pixel_width(0.0001, 1000.0);
  ray::OptionalKeyValueArgument pixel_width_option("pixel_width", 'p', &pixel_width);
  ray::FileArgument cloud_file, image_file;
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &viewpoint, &style, &image_file}, {&pixel_width_option}))
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();

  double pix_width = pixel_width_option.isSet() ? pixel_width.value() : cloud.estimatePointSpacing();
  Eigen::Vector3d min_bounds = cloud.calcMinBound();
  Eigen::Vector3d max_bounds = cloud.calcMaxBound();
  Eigen::Vector3d extent = max_bounds - min_bounds;
  int axis = 0;
  if (viewpoint.selectedKey() == "top")
    axis = 2;
  else if (viewpoint.selectedKey() == "front" || viewpoint.selectedKey() == "back")
    axis = 1;
  double dir = 1;
  if (viewpoint.selectedKey() == "left" || viewpoint.selectedKey() == "front")
    dir = -1;
  
  int ax1 = (axis+1)%3;
  int ax2 = (axis+2)%3;
  int width  = 1 + static_cast<int>(extent[ax1] / pix_width);
  int height = 1 + static_cast<int>(extent[ax2] / pix_width);
  std::cout << "outputting " << width << "x" << height << " image" << std::endl;
  std::vector<Eigen::Vector4d> pixels(width * height); 
  memset(&pixels[0], 0, sizeof(Eigen::Vector4d) * width*height);
  double max_val = 0.0;
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    ray::RGBA &colour = cloud.colours[i];
    Eigen::Vector3d col(colour.red, colour.green, colour.blue);
    Eigen::Vector3d point = style.selectedID() == 4 ? cloud.starts[i] : cloud.ends[i];
    Eigen::Vector3d pos = (point - min_bounds) / pix_width;
    Eigen::Vector3i p = (pos).cast<int>();
    int x = p[ax1], y = p[ax2];
    Eigen::Vector4d &pix = pixels[x + width*y];
    switch (style.selectedID())
    {
      case 0: // ends
      case 4: // starts
        if (pos[axis]*dir > pix[3]*dir || pix[3] == 0.0)
          pix = Eigen::Vector4d(col[0], col[1], col[2], pos[axis]);
        break;
      case 1: // mean
        pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
        break;
      case 2: // sum
        pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
        max_val = std::max(max_val, pix[0]);
        max_val = std::max(max_val, pix[1]);
        max_val = std::max(max_val, pix[2]);
        break;
      case 3: // density
        break;
      case 5: // rays
        // walk each ray:
        break;
      default:
        break;
    }
  }
  std::vector<ray::RGBA> pixel_colours;
  std::vector<Eigen::Vector3f> float_pixel_colours;
  bool is_hdr = image_file.nameExt() == ".hdr";
  if (is_hdr)
    float_pixel_colours.resize(width * height);
  else
    pixel_colours.resize(width*height);

  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Eigen::Vector4d colour = pixels[x + width*y];
      Eigen::Vector3d col3d(colour[0], colour[1], colour[2]);
      switch (style.selectedID())
      {
        case 1: // mean
          col3d /= colour[3];
          break;
        case 2: // sum
          col3d *= 255.0 / max_val;
          break;
        default:
          break;
      }
      if (is_hdr)
        float_pixel_colours[x + width * y] = col3d.cast<float>();
      else
      {
        ray::RGBA col;
        col.red   = uint8_t(std::min(col3d[0], 255.0));
        col.green = uint8_t(std::min(col3d[1], 255.0));
        col.blue  = uint8_t(std::min(col3d[2], 255.0));
        col.alpha = 255;
        pixel_colours[x + width * (height - 1 - y)] = col;
      }
    }
  }
  const char *image_name = image_file.name().c_str();
  if (image_file.nameExt() == "png")
    stbi_write_png(image_name, width, height, 4, (void *)&pixel_colours[0], 4 * width);
  else if (image_file.nameExt() == "bmp")
    stbi_write_bmp(image_name, width, height, 4, (void *)&pixel_colours[0]);
  else if (image_file.nameExt() == "tga")
    stbi_write_tga(image_name, width, height, 4, (void *)&pixel_colours[0]);
  else if (image_file.nameExt() == "png")
    stbi_write_jpg(image_name, width, height, 4, (void *)&pixel_colours[0], 100); // maximal quality
  else if (image_file.nameExt() == "hdr")
    stbi_write_hdr(image_name, width, height, 4, (float *)&float_pixel_colours[0]);
  else
  {
    std::cerr << "Error: image format " << image_file.nameExt() << " not known" << std::endl;
    usage();
  }

  return 0;
}
