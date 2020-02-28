#pragma once
#include "rayutils.h"
#include "raypose.h"
namespace RAY
{

struct Trajectory
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

    /// Operators
    Node operator *(double scale) const
    {
      return Node(pose*scale, time*scale);
    }
    Node operator /(double scale) const
    {
      ASSERT(scale != 0.0);
      return Node(pose / scale, time / scale);
    }
    Node operator *(const Pose &otherPose) const
    {
      return Node(pose*otherPose, time);
    }
    Node operator +(const Node &node) const
    {
      return Node(pose + node.pose, time + node.time);
    }
    Node operator -(const Node &node) const
    {
      return Node(pose - node.pose, time - node.time);
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