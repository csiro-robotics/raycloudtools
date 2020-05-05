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

namespace RAY
{

struct RAYLIB_EXPORT Trajectory
{
  struct Node
  {
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

  void save(const std::string &fileName, double timeOffset = 0.0);
  bool load(const std::string &fileName);
};


}

#endif // RAYLIB_RAYTRAJECTORY_H
