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

class RAYLIB_EXPORT Trajectory
{
public:
  class Node
  {
  public:
    Pose pose;
    double time;

    Node() {}
    Node(const Pose &pose, double time)
    {
      this->pose = pose;
      this->time = time;
    }
    
    inline void normalise()
    {
      pose.normalise();
    }
  };
  std::vector<Node> nodes;

  void save(const std::string &file_name, double time_offset = 0.0);
  bool load(const std::string &file_name);
};


}

#endif // RAYLIB_RAYTRAJECTORY_H
