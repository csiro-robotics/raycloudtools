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
/// @brief subsample to 1 point per @c vox_width wide voxel in metres
/// This is a spatially even subsampling, but also emphasises outlier as a side-effect
bool RAYLIB_EXPORT decimateSpatial(const std::string &file_stub, double vox_width);

/// @brief subsample to every @c num_rays rays
/// This is an unbiased subsampling, but will be over-sampled in stationary areas as a side-effect
/// Note that while this is called temporal decimation, it decimates evenly in file order, which isn't 
/// necessarily temporal order. Though it typically is stored that way on single scans.
bool RAYLIB_EXPORT decimateTemporal(const std::string &file_stub, int num_rays);

/// @brief subsample to @c num_rays rays (temporally decimated) for each @c vox_width wide voxel
/// This allows a more even distribution of points while maintaining details better than pure spatial decimation
bool RAYLIB_EXPORT decimateSpatioTemporal(const std::string &file_stub, double vox_width, int num_rays);

/// @brief decimate to no more than 1 point per voxel of width @c radius_per_length x ray length. 
/// This is used when error is proportional to ray length, prioritising closer measurements and leaving distant areas sparse
bool RAYLIB_EXPORT decimateAngular(const std::string &file_stub, double radius_per_length);
}  // namespace ray

#endif  // RAYLIB_RAYDECIMATION_H
