// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTREEGEN_H
#define RAYLIB_RAYTREEGEN_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"

namespace RAY
{
#define pitchAngle (30.0 * pi/180.0) // alternative to above, defined by relative angle of the head of the two branches
#define splitAngle (45.0 * pi/180.0) // angle in the Y shape, usually around 45 degrees
#define branchGradient 20.0 // length per radius

double RAYLIB_EXPORT getMainBranchAngle(double covarianceAngle);
void RAYLIB_EXPORT fillBranchAngleLookup();

struct RAYLIB_EXPORT TreeGen
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

#endif // RAYLIB_RAYTREEGEN_H
