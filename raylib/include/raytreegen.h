#pragma once
#include "rayutils.h"

namespace RAY
{
struct TreeGen
{
  vector<Vector3d> leaves;
  vector<Vector3d> rayStarts, rayEnds;
  struct Branch
  {
    Vector3d tip;
    double radius;
    int parentIndex;
  };
  vector<Branch> branches;
  Vector3d root;

  void addBranch(int parentIndex, Pose pose, double radius, double randomFactor);
  // create the tree structure, and list of leaf points
  void make(const Vector3d &rootPos, double trunkRadius, double randomFactor = 0.0); // 0 to 1
  // create a set of rays covering the tree at a roughly uniform distribution
  void generateRays(double rayDensity);
};

}
