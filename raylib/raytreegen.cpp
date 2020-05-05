// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytreegen.h"
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

using namespace RAY;
using namespace std;
using namespace Eigen;

void getBranchInfo(double mainBranchAngle, double &secondaryBranchAngle, double &mainBranchRadius, double &secondaryBranchRadius)
{
  secondaryBranchAngle = splitAngle - mainBranchAngle;
  double sinang1 = sin(mainBranchAngle);
  double sinang2 = sin(secondaryBranchAngle);
  double q1 = sinang1/sinang2;
  double q2 = sinang2/sinang1;
  mainBranchRadius = 1.0 / sqrt(sqrt(q1)+1.0);
  secondaryBranchRadius = 1.0 / sqrt(sqrt(q2)+1.0);
}

static const int branchAngleSize = 400;
static double branchAngleLookup[branchAngleSize];
double RAY::getMainBranchAngle(double covarianceAngle)
{
  double x = covarianceAngle * (double)branchAngleSize / (0.5*pi);
  int i = max(0, min((int)x, branchAngleSize-2));
  double blend = max(0.0, min(x - (double)i, 1.0));
  return branchAngleLookup[i]*(1.0-blend) + branchAngleLookup[i+1]*blend;
}

void RAY::fillBranchAngleLookup()
{
  double lastG = 0.0;
  double lastAng1 = 0.0;
  int j = 0;
  branchAngleLookup[j++] = 0.0;
  double gStep = (pi/2.0) / (double)branchAngleSize;
 // int numSamples = 100*branchAngleSize;
  double ang1 = 4e-10;
  double scale = 2.0;
  while (ang1 != splitAngle*0.5)
  {
    ang1 *= 1.0 + scale;
    ang1 = min(ang1, splitAngle*0.5);
  //  double ang1 = splitAngle*0.5 * (double)i/(double)(numSamples-1);
    double ang2, l1, l2;
    getBranchInfo(ang1, ang2, l1, l2);
    double g = atan2(l1*sin(ang1) + l2*sin(ang2), l1*cos(ang1) - l2*cos(ang2));
    // i.e. we have g for each ang1e of the branch... we want to go the other way....
    int nums = 0;
    while (g > (double)j*gStep && j<branchAngleSize)
    {
      double blend = ((double)j*gStep - lastG) / (g-lastG);
      branchAngleLookup[j++] = lastAng1*(1.0 - blend) + ang1*blend;
      nums++;
    }
    if (nums == 0)
      scale *= 1.5;
    if (nums > 1)
      scale /= (double)nums;
    lastG = g;
    lastAng1 = ang1;
  }
}

static const double minimumRadius = 0.001;
static Vector3d COM(0,0,0);
static double totalMass = 0.0;

void TreeGen::addBranch(int parentIndex, Pose pose, double radius, double randomFactor)
{
  if (radius < minimumRadius)
  {
    leaves.push_back(pose.position);
    return;
  }
  Vector3d p1 = pose.position;
  double randScale = random(1.0-randomFactor, 1.0+randomFactor);
  pose.position = pose * Vector3d(0, 0, radius*branchGradient*randScale);
  double phi = (sqrt(5) + 1.0)/2.0;
  randScale = random(1.0-randomFactor, 1.0+randomFactor);
  pose.rotation = pose.rotation * Quaterniond(AngleAxisd(2.0*pi*phi*randScale, Vector3d(0,0,1)));
  Branch branch;
  branch.tip = pose.position;
  branch.radius = radius;
  branch.parentIndex = parentIndex;
  int index = int(branches.size());
  branches.push_back(branch);
  COM += radius*radius*radius*(branch.tip + p1)/2.0;
  totalMass += radius*radius*radius;
  
  Pose child1 = pose, child2 = pose;
  double angle1, angle2, radius1, radius2;
  randScale = random(1.0-randomFactor, 1.0+randomFactor);
  angle1 = getMainBranchAngle(pi*0.5 - pitchAngle*randScale); // splitAngle*0.5 * mainBranchAngleRatio;
  getBranchInfo(angle1, angle2, radius1, radius2);
  randScale = random(1.0-randomFactor, 1.0+randomFactor);
  angle2 *= randScale;
  radius1 *= radius;
  radius2 *= radius;

  Quaterniond q1(AngleAxisd(angle1, Vector3d(1,0,0)));
  Quaterniond q2(AngleAxisd(angle2, Vector3d(-1,0,0)));
  
  child1.rotation = child1.rotation * q1;
  child2.rotation = child2.rotation * q2;
  addBranch(index, child1, radius1, randomFactor);
  addBranch(index, child2, radius2, randomFactor);
}

// create the tree structure, and list of leaf points
void TreeGen::make(const Vector3d &rootPos, double trunkRadius, double randomFactor)
{
  COM.setZero();
  totalMass = 0.0;

  Pose base(rootPos, Quaterniond(AngleAxisd(randomFactor*random(0.0, 2.0*pi), Vector3d(0, 0, 1))));
  Branch branch;
  branch.tip = rootPos;
  branch.parentIndex = -1;
  branch.radius = 2.0*trunkRadius;
  branches.push_back(branch);
  addBranch(0, base, trunkRadius, randomFactor);

  COM /= totalMass;
 // cout << "COM: " << COM.transpose() << ", grad = " << COM[2]/trunkRadius << endl;
  double scale = branchGradient/(COM[2]/trunkRadius);
  for (auto &leaf: leaves)
    leaf = rootPos + (leaf-rootPos)*scale;
  for (auto &start: rayStarts)
    start = rootPos + (start-rootPos)*scale;
  for (auto &end: rayEnds)
    end = rootPos + (end - rootPos)*scale;
  for (auto &branch: branches)
    branch.tip = rootPos + (branch.tip - rootPos)*scale;
}

// create a set of rays covering the tree at a roughly uniform distribution
void TreeGen::generateRays(double rayDensity)
{
  vector<double> cumulativeSize(branches.size());
  cumulativeSize[0] = 0;
  for (int i = 1; i<(int)branches.size(); i++)
  {
    Branch &branch = branches[i];
    Branch &parentBranch = branches[branch.parentIndex];
    double area = (branch.tip - parentBranch.tip).norm() * 2.0*pi*(branch.radius + parentBranch.radius)/2.0; // slightly approximate
    cumulativeSize[i] = cumulativeSize[i-1] + area;
  }
  
  int numRays = (int)(rayDensity * cumulativeSize.back());
  double areaPerRay = cumulativeSize.back() / (double)numRays;
  double totalArea = 0.0;
  int i = 0;
  for (int j = 0; j<numRays; j++)
  {
    totalArea += areaPerRay;
    while (i<(int)cumulativeSize.size()-1 && cumulativeSize[i] < totalArea)
      i++;
  
    Branch &branch = branches[i];
    Branch &parentBranch = branches[branch.parentIndex];
    
    // simplest is to randomise a point on a cone.... it will have overlap and gaps, but it is a starting point...
    double r1 = branch.radius;
    double r2 = branch.radius + 0.1*(parentBranch.radius - branch.radius);            
    double r = sqrt(random(sqr(r1), sqr(r2)));
    double t = (r - r1)/(r2-r1);
    Vector3d online = branch.tip + (parentBranch.tip - branch.tip)*t;
    double angle = random(0.0, 2.0*pi);
    Vector3d up = (branch.tip - parentBranch.tip).normalized();
    Vector3d side = up.cross(Vector3d(1,2,3));
    Vector3d fwd = up.cross(side).normalized();
    side = fwd.cross(up);
    Vector3d offset = side * sin(angle) + fwd*cos(angle);
    Vector3d pos = online + offset * r;
    if (!(pos[0] == pos[0]))
      cout << "bad pos" << pos.transpose() << endl;
    rayEnds.push_back(pos);
    Vector3d from = Vector3d(random(-1,1), random(-1,1), random(-1,1));
    if (from.dot(offset) < 0.0)
      from = -from;
    rayStarts.push_back(pos + from*0.1*r1);
  }
}
