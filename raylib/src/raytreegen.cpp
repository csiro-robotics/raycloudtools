/* (c) Copyright CSIRO 2013. Author: Thomas Lowe
   This software is provided under the terms of Schedule 1 of the license agreement between CSIRO, 3DLM and GeoSLAM.
*/
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include "raytreegen.h"
using namespace Ray;
using namespace std;
using namespace Eigen;

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
  pose.rotation = pose.rotation * Quat(Vector3d(0,0,2.0*pi*phi*randScale));
  Branch branch;
  branch.tip = pose.position;
  branch.radius = radius;
  branch.parentIndex = parentIndex;
  int index = branches.size();
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

  Quat q1(Vector3d(angle1, 0,0));
  Quat q2(Vector3d(-angle2,0,0));
  
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

  Pose base(rootPos, Quat(Vector3d(0, 0, randomFactor*random(0.0, 2.0*pi))));
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
