// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYALIGNMENT_H
#define RAYLIB_RAYALIGNMENT_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"

#include <complex>

typedef std::complex<double> Complex;

namespace ray
{
/// Coarse raycloud alignment. This translates and 'yaw's the ray cloud, under the common assumption that pitch and roll
/// are already accurate. Transforms the first cloud in the pair @c clouds, to align with the second.
/// This is a cross-correlation method that requires a @c voxel_width (typically on the order of a metre)
/// the @c verbose argument saves out plan-view images at each step of the method.
/// The method uses a scale-free Fourier-Mellin transform to efficiently cross-correlate the cloud's end point densities.
void RAYLIB_EXPORT alignCloud0ToCloud1(Cloud *clouds, double voxel_width, bool verbose = false);
}  // namespace ray

#endif  // RAYLIB_RAYALIGNMENT_H
