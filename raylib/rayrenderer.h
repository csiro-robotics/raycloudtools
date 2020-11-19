// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTRAJECTORY_H
#define RAYLIB_RAYTRAJECTORY_H

#include "raylib/raylibconfig.h"
#include "raycuboid.h"
#include "rayutils.h"
#include "raypose.h"

namespace ray
{
/// Supported view directions on cloud data
enum class RAYLIB_EXPORT ViewDirection
{
  Top, 
  Left, 
  Right, 
  Front,
  Back
};

/// Supported render styles on cloud data
enum class RAYLIB_EXPORT RenderStyle
{
  Ends, 
  Mean, 
  Sum, 
  Starts, 
  Rays, 
  Density, 
  Density_rgb
};

/// Render a ray cloud according to the supplies parameters
bool RAYLIB_EXPORT renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction, 
                               RenderStyle style, double pix_width, const std::string &image_file);
}
#endif // RAYLIB_RAYTRAJECTORY_H
