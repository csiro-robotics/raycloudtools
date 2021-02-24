// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTERRAIN_H
#define RAYLIB_RAYTERRAIN_H

#include "raylib/raylibconfig.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include "../rayalignment.h"
typedef Eigen::Matrix<double, 4, 1> Vector4d;

namespace ray
{
  
/// A class for storage and extraction of a smooth ground height function
class RAYLIB_EXPORT Terrain
{
public:
  /// Extracts a robust smooth surface from the ray cloud, into the state
  void extract(const Cloud &cloud, const std::string &file_prefix, double gradient, bool verbose);
private:
  void getParetoFront(const std::vector<Vector4d> &points, std::vector<Vector4d> &front);
};


}  // namespace ray

#endif  // RAYLIB_RAYTERRAIN_H
