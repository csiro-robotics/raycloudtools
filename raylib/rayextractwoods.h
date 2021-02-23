// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYEXTRACT_WOODS_H
#define RAYLIB_RAYEXTRACT_WOODS_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"


namespace ray
{
#define LENGTH

struct Ray
{
  Ray(const Eigen::Vector3d &start, const Eigen::Vector3d &pos) : start(start), pos(pos) {}
  Eigen::Vector3d start, pos;
};

struct Cell
{
  std::vector<Ray> rays;
  Eigen::Vector2d minBound;
  double height;
};

struct Wood
{
  Wood(const Cloud &cloud, double midRadius, bool verbose);
  bool save(const std::string &filename);
  std::vector<Trunk> trunk_bases;
};

struct Accumulator
{
  Accumulator(): weight(0), x(0), abs_x(0), y(0,0), xy(0,0), x2(0), radius(0), radius2(0), z(0,0), xz(0,0) {}

  double weight;
  double x;
  double abs_x;
  Eigen::Vector2d y;
  Eigen::Vector2d xy;
  double x2;
  double radius;
  double radius2;
  Eigen::Vector2d z;
  Eigen::Vector2d xz;
};
} // namespace ray
#endif // RAYLIB_RAYEXTRACT_WOODS_H