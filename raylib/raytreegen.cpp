// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytreegen.h"
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

using namespace ray;
using namespace std;
using namespace Eigen;

void getBranchInfo(double main_branch_angle, double &secondary_branch_angle, double &main_branch_radius, double &secondary_branch_radius)
{
  secondary_branch_angle = splitAngle - main_branch_angle;
  double sinang1 = sin(main_branch_angle);
  double sinang2 = sin(secondary_branch_angle);
  double q1 = sinang1/sinang2;
  double q2 = sinang2/sinang1;
  main_branch_radius = 1.0 / sqrt(sqrt(q1)+1.0);
  secondary_branch_radius = 1.0 / sqrt(sqrt(q2)+1.0);
}

static const int kBranchAngleSize = 400;
static double branch_angle_lookup[kBranchAngleSize];
double ray::getMainBranchAngle(double covariance_angle)
{
  double x = covariance_angle * (double)kBranchAngleSize / (0.5*kPi);
  int i = max(0, min((int)x, kBranchAngleSize-2));
  double blend = max(0.0, min(x - (double)i, 1.0));
  return branch_angle_lookup[i]*(1.0-blend) + branch_angle_lookup[i+1]*blend;
}

void ray::fillBranchAngleLookup()
{
  double last_g = 0.0;
  double last_ang1 = 0.0;
  int j = 0;
  branch_angle_lookup[j++] = 0.0;
  double g_step = (kPi/2.0) / (double)kBranchAngleSize;
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
    while (g > (double)j*g_step && j<kBranchAngleSize)
    {
      double blend = ((double)j*g_step - last_g) / (g-last_g);
      branch_angle_lookup[j++] = last_ang1*(1.0 - blend) + ang1*blend;
      nums++;
    }
    if (nums == 0)
      scale *= 1.5;
    if (nums > 1)
      scale /= (double)nums;
    last_g = g;
    last_ang1 = ang1;
  }
}

static const double kMinimumRadius = 0.001;
static Vector3d com(0,0,0);
static double total_mass = 0.0;

void TreeGen::addBranch(int parent_index, Pose pose, double radius, double random_factor)
{
  if (radius < kMinimumRadius)
  {
    leaves.push_back(pose.position);
    return;
  }
  Vector3d p1 = pose.position;
  double rand_scale = random(1.0-random_factor, 1.0+random_factor);
  pose.position = pose * Vector3d(0, 0, radius*branchGradient*rand_scale);
  double phi = (sqrt(5) + 1.0)/2.0;
  rand_scale = random(1.0-random_factor, 1.0+random_factor);
  pose.rotation = pose.rotation * Quaterniond(AngleAxisd(2.0*kPi*phi*rand_scale, Vector3d(0,0,1)));
  Branch branch;
  branch.tip = pose.position;
  branch.radius = radius;
  branch.parent_index = parent_index;
  int index = int(branches.size());
  branches.push_back(branch);
  com += radius*radius*radius*(branch.tip + p1)/2.0;
  total_mass += radius*radius*radius;
  
  Pose child1 = pose, child2 = pose;
  double angle1, angle2, radius1, radius2;
  rand_scale = random(1.0-random_factor, 1.0+random_factor);
  angle1 = getMainBranchAngle(kPi*0.5 - pitchAngle*rand_scale); // splitAngle*0.5 * mainBranchAngleRatio;
  getBranchInfo(angle1, angle2, radius1, radius2);
  rand_scale = random(1.0-random_factor, 1.0+random_factor);
  angle2 *= rand_scale;
  radius1 *= radius;
  radius2 *= radius;

  Quaterniond q1(AngleAxisd(angle1, Vector3d(1,0,0)));
  Quaterniond q2(AngleAxisd(angle2, Vector3d(-1,0,0)));
  
  child1.rotation = child1.rotation * q1;
  child2.rotation = child2.rotation * q2;
  addBranch(index, child1, radius1, random_factor);
  addBranch(index, child2, radius2, random_factor);
}

// create the tree structure, and list of leaf points
void TreeGen::make(const Vector3d &root_pos, double trunk_radius, double random_factor)
{
  com.setZero();
  total_mass = 0.0;

  Pose base(root_pos, Quaterniond(AngleAxisd(random_factor*random(0.0, 2.0*kPi), Vector3d(0, 0, 1))));
  Branch branch;
  branch.tip = root_pos;
  branch.parent_index = -1;
  branch.radius = 2.0*trunk_radius;
  branches.push_back(branch);
  addBranch(0, base, trunk_radius, random_factor);

  com /= total_mass;
 // cout << "COM: " << COM.transpose() << ", grad = " << COM[2]/trunkRadius << endl;
  double scale = branchGradient/(com[2]/trunk_radius);
  for (auto &leaf: leaves)
    leaf = root_pos + (leaf-root_pos)*scale;
  for (auto &start: ray_starts)
    start = root_pos + (start-root_pos)*scale;
  for (auto &end: ray_ends)
    end = root_pos + (end - root_pos)*scale;
  for (auto &branch: branches)
    branch.tip = root_pos + (branch.tip - root_pos)*scale;
}

// create a set of rays covering the tree at a roughly uniform distribution
void TreeGen::generateRays(double ray_density)
{
  vector<double> cumulative_size(branches.size());
  cumulative_size[0] = 0;
  for (int i = 1; i<(int)branches.size(); i++)
  {
    Branch &branch = branches[i];
    Branch &parent_branch = branches[branch.parent_index];
    double area = (branch.tip - parent_branch.tip).norm() * 2.0*kPi*(branch.radius + parent_branch.radius)/2.0; // slightly approximate
    cumulative_size[i] = cumulative_size[i-1] + area;
  }
  
  int num_rays = (int)(ray_density * cumulative_size.back());
  double area_per_ray = cumulative_size.back() / (double)num_rays;
  double total_area = 0.0;
  int i = 0;
  for (int j = 0; j<num_rays; j++)
  {
    total_area += area_per_ray;
    while (i<(int)cumulative_size.size()-1 && cumulative_size[i] < total_area)
      i++;
  
    Branch &branch = branches[i];
    Branch &parent_branch = branches[branch.parent_index];
    
    // simplest is to randomise a point on a cone.... it will have overlap and gaps, but it is a starting point...
    double r1 = branch.radius;
    double r2 = branch.radius + 0.1*(parent_branch.radius - branch.radius);            
    double r = sqrt(random(sqr(r1), sqr(r2)));
    double t = (r - r1)/(r2-r1);
    Vector3d online = branch.tip + (parent_branch.tip - branch.tip)*t;
    double angle = random(0.0, 2.0*kPi);
    Vector3d up = (branch.tip - parent_branch.tip).normalized();
    Vector3d side = up.cross(Vector3d(1,2,3));
    Vector3d fwd = up.cross(side).normalized();
    side = fwd.cross(up);
    Vector3d offset = side * sin(angle) + fwd*cos(angle);
    Vector3d pos = online + offset * r;
    if (!(pos[0] == pos[0]))
      cout << "bad pos" << pos.transpose() << endl;
    ray_ends.push_back(pos);
    Vector3d from = Vector3d(random(-1,1), random(-1,1), random(-1,1));
    if (from.dot(offset) < 0.0)
      from = -from;
    ray_starts.push_back(pos + from*0.1*r1);
  }
}
