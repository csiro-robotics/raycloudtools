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
  void extract(const Cloud &cloud, const std::string &file_prefix, bool verbose);
private:
  inline bool greaterThan(const Eigen::Vector4d &a, const Eigen::Vector4d &b, bool positive)
  {
    if (positive)
      return a[0] >= b[0] && a[1] >= b[1] && a[2] >= b[2];
    return a[0] <= b[0] && a[1] <= b[1] && a[2] <= b[2];
  }
  void getParetoFront(std::vector<Vector4d> points, std::vector<Vector4d> &front, bool pos);
  void getParetoFrontFast(std::vector<Vector4d> points, std::vector<Vector4d> &front, bool pos);
};


}  // namespace ray

#endif  // RAYLIB_RAYTERRAIN_H
