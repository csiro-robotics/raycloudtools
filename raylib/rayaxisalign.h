// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYAXISALIGN_H
#define RAYLIB_RAYAXISALIGN_H

#include <string>
#include "raylib/raylibconfig.h"

namespace ray
{
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

#endif  // RAYLIB_RAYAXISALIGN_H
