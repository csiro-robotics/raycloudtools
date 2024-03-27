// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tim Devereux
#ifndef RAYLIB_raydensegrid_H
#define RAYLIB_raydensegrid_H

#include "raylib/raylibconfig.h"
#include "raylib/rayutils.h"

namespace ray
{
bool generateAreaVoxels(const std::string &cloud_stub, double voxel_size);
}  // namespace ray
#endif  // RAYLIB_RAYLEAVES_H
