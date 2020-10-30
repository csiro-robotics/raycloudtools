// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTRAJECTORY_H
#define RAYLIB_RAYTRAJECTORY_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"

namespace ray
{

/// Simple class to represent a moving point over time
class RAYLIB_EXPORT Trajectory
{
public:
  class RAYLIB_EXPORT Node
  {
  public:
    Eigen::Vector3d pos;
    double time;

    Node() {}
    Node(const Eigen::Vector3d &pos, double time)
    {
      this->pos = pos;
      this->time = time;
    }
  };
  std::vector<Node> nodes;

  /// Save trajectory to a text file. One line per Node
  void save(const std::string &file_name);

  /// Load trajectory from file. The file is expected to be a text file, with one Node entry per line
  bool load(const std::string &file_name);

  void calculateStartPoints(const std::vector<double> &times, std::vector<Eigen::Vector3d> &starts);
};


}

#endif // RAYLIB_RAYTRAJECTORY_H
