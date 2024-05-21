// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYDECIMATION_H
#define RAYLIB_RAYDECIMATION_H

#include <iostream>
#include <limits>
#include "raycloud.h"
#include "rayutils.h"

namespace ray
{
bool RAYLIB_EXPORT decimateSpatial(const std::string &file_stub, double vox_width);
bool RAYLIB_EXPORT decimateTemporal(const std::string &file_stub, int num_rays);
bool RAYLIB_EXPORT decimateSpatioTemporal(const std::string &file_stub, double vox_width, int num_rays);
bool RAYLIB_EXPORT decimateAngular(const std::string &file_stub, double radius_per_length);
}  // namespace ray

#endif  // RAYLIB_RAYDECIMATION_H
