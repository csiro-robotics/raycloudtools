// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYSPLITTER_H
#define RAYLIB_RAYSPLITTER_H

#include <iostream>
#include <limits>
#include "raycloud.h"
#include "rayutils.h"

namespace ray
{
/// Split a file into @c in_name or @c out_name depending on the function @c is_outside.
bool RAYLIB_EXPORT split(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                         std::function<bool(const Cloud &cloud, int i)> is_outside);

/// Split a ray cloud around a plane. This splits individual rays, to maintain the validity of the cloud
/// Each ray goes into file @c in_name or @c out_name depending on the side of the plane
bool RAYLIB_EXPORT splitPlane(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                              const Eigen::Vector3d &plane);

/// Split a ray cloud around a cuboid defined by @c centre and @c extents. This also splits individual rays.
/// The results go into file @c in_name or @c out_name depending on which side of the box each ray is on
/// With @c in_name becoming the cloud cropped to the bounding box
bool RAYLIB_EXPORT splitBox(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                            const Eigen::Vector3d &centre, const Eigen::Vector3d &extents);

/// Split a ray cloud into a grid of files, named with suffix _X_Y_Z.ply, for each grid coordinate X,Y,Z.
/// Aligned so that cell 0,0,0 is centred at 0,0,0, and has dimensions @c cell_width
/// @c overlap generates larger cells so that they overlap by the specified value
bool RAYLIB_EXPORT splitGrid(const std::string &file_name, const std::string &cloud_name_stub,
                             const Eigen::Vector3d &cell_width, double overlap = 0.0);

/// Split a ray cloud into a grid of files, named with suffix _X_Y_Z_T.ply, for each grid coordinate X,Y,Z,T.
/// Aligned so that cell 0,0,0,0 is centred at 0,0,0,0 and has dimensions @c cell_width
/// @c overlap generates larger cells so that they overlap by the specified value
bool RAYLIB_EXPORT splitGrid(const std::string &file_name, const std::string &cloud_name_stub,
                             const Eigen::Vector4d &cell_width, double overlap = 0.0);

/// Split a ray cloud into one cloud per colour, ignoring differences in alpha. For example, when identified objects in
/// the cloud are given a unique colour
bool RAYLIB_EXPORT splitColour(const std::string &file_name, const std::string &cloud_name_stub);

}  // namespace ray

#endif  // RAYLIB_RAYSPLITTER_H
