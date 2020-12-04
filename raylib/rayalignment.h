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
/// NOTE @c clouds is a pair of clouds, it should point to an array with at least 2 elements
void RAYLIB_EXPORT alignCloud0ToCloud1(Cloud *clouds, double voxel_width, bool verbose = false);

/// Align cloud to its principle axes using 3D translation and 2D rotation (yaw about Z). 
/// The highest density orthogonal planes define the new coordinate frame.
/// The cloud is therefore translated to the intersection of these planes and 'yawed' such that the strongest 
/// vertical plane represents the y axis. 
/// Ambiguity with respect to 180 degrees is resolved as follows:
/// if the centroid y component is farther from the origin than its x component, then make centroid y positive
/// otherwise make centroid x positive.
/// The purpose of this is to make the chosen alignment as robust as possible to variation, or rescans.
bool RAYLIB_EXPORT alignCloudToAxes(const std::string &cloud_name, const std::string &aligned_file);
}  // namespace ray

#endif  // RAYLIB_RAYALIGNMENT_H
