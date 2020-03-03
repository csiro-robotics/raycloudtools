#pragma once
#include "rayutils.h"
#include "raypose.h"

namespace RAY
{
struct TreeGen
{
  std::vector<Eigen::Vector3d> leaves;
  std::vector<Eigen::Vector3d> rayStarts, rayEnds;
  struct Branch
  {
    Eigen::Vector3d tip;
    double radius;
    int parentIndex;
  };
  std::vector<Branch> branches;
  Eigen::Vector3d root;

  void addBranch(int parentIndex, Pose pose, double radius, double randomFactor);
  // create the tree structure, and list of leaf points
  void make(const Eigen::Vector3d &rootPos, double trunkRadius, double randomFactor = 0.0); // 0 to 1
  // create a set of rays covering the tree at a roughly uniform distribution
  void generateRays(double rayDensity);
};

}
