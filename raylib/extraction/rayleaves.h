// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLEAVES_H
#define RAYLIB_RAYLEAVES_H

#include "raylib/raylibconfig.h"
#include "raylib/rayutils.h"

namespace ray
{
bool generateLeaves(const std::string &cloud_stub, const std::string &trees_file, const std::string &leaf_file,
                    double leaf_area, double droop, int distribution, double leafAreaDensity, bool stalks);
}  // namespace ray
#endif  // RAYLIB_RAYLEAVES_H
