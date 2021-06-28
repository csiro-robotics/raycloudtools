// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytreegen.h"
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

namespace ray
{
void getBranchInfo(double main_branch_angle, double &secondary_branch_angle, double &main_branch_radius,
                   double &secondary_branch_radius)
{
  secondary_branch_angle = splitAngle - main_branch_angle;
  double sinang1 = sin(main_branch_angle);
  double sinang2 = sin(secondary_branch_angle);
  double q1 = sinang1 / sinang2;
  double q2 = sinang2 / sinang1;
  main_branch_radius = 1.0 / sqrt(sqrt(q1) + 1.0);
  secondary_branch_radius = 1.0 / sqrt(sqrt(q2) + 1.0);
}

static const int kBranchAngleSize = 400;
static double branch_angle_lookup[kBranchAngleSize];
double getMainBranchAngle(double covariance_angle)
{
  double x = covariance_angle * (double)kBranchAngleSize / (0.5 * kPi);
  int i = std::max(0, std::min((int)x, kBranchAngleSize - 2));
  double blend = std::max(0.0, std::min(x - (double)i, 1.0));
  return branch_angle_lookup[i] * (1.0 - blend) + branch_angle_lookup[i + 1] * blend;
}

void fillBranchAngleLookup()
{
  double last_g = 0.0;
  double last_ang1 = 0.0;
  int j = 0;
  branch_angle_lookup[j++] = 0.0;
  double g_step = (kPi / 2.0) / (double)kBranchAngleSize;
  // int numSamples = 100*branchAngleSize;
  double ang1 = 4e-10;
  double scale = 2.0;
  while (ang1 != splitAngle * 0.5)
  {
    ang1 *= 1.0 + scale;
    ang1 = std::min(ang1, splitAngle * 0.5);
    //  double ang1 = splitAngle*0.5 * (double)i/(double)(numSamples-1);
    double ang2, l1, l2;
    getBranchInfo(ang1, ang2, l1, l2);
    double g = atan2(l1 * sin(ang1) + l2 * sin(ang2), l1 * cos(ang1) - l2 * cos(ang2));
    // i.e. we have g for each ang1e of the branch... we want to go the other way....
    int nums = 0;
    while (g > (double)j * g_step && j < kBranchAngleSize)
    {
      double blend = ((double)j * g_step - last_g) / (g - last_g);
      branch_angle_lookup[j++] = last_ang1 * (1.0 - blend) + ang1 * blend;
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
static Eigen::Vector3d com(0, 0, 0);
static double total_mass = 0.0;

void TreeGen::addBranch(int parent_index, Pose pose, double radius, double random_factor)
{
  if (radius < kMinimumRadius)
  {
    leaves_.push_back(pose.position);
    return;
  }
  Eigen::Vector3d p1 = pose.position;
  double rand_scale = random(1.0 - random_factor, 1.0 + random_factor);
  pose.position = pose * Eigen::Vector3d(0, 0, radius * branchGradient * rand_scale);
  double phi = (sqrt(5) + 1.0) / 2.0;
  rand_scale = random(1.0 - random_factor, 1.0 + random_factor);
  pose.rotation = pose.rotation * Eigen::Quaterniond(Eigen::AngleAxisd(2.0 * kPi * phi * rand_scale, Eigen::Vector3d(0, 0, 1)));
  Branch branch;
  branch.tip = pose.position;
  branch.radius = radius;
  branch.parent_index = parent_index;
  int index = int(branches_.size());
  branches_.push_back(branch);
  com += radius * radius * radius * (branch.tip + p1) / 2.0;
  total_mass += radius * radius * radius;

  Pose child1 = pose, child2 = pose;
  double angle1, angle2, radius1, radius2;
  rand_scale = random(1.0 - random_factor, 1.0 + random_factor);
  angle1 = getMainBranchAngle(kPi * 0.5 - pitchAngle * rand_scale);  // splitAngle*0.5 * mainBranchAngleRatio;
  getBranchInfo(angle1, angle2, radius1, radius2);
  rand_scale = random(1.0 - random_factor, 1.0 + random_factor);
  angle2 *= rand_scale;
  radius1 *= radius;
  radius2 *= radius;

  Eigen::Quaterniond q1(Eigen::AngleAxisd(angle1, Eigen::Vector3d(1, 0, 0)));
  Eigen::Quaterniond q2(Eigen::AngleAxisd(angle2, Eigen::Vector3d(-1, 0, 0)));

  child1.rotation = child1.rotation * q1;
  child2.rotation = child2.rotation * q2;
  addBranch(index, child1, radius1, random_factor);
  addBranch(index, child2, radius2, random_factor);
}

// create the tree structure, and list of leaf points
void TreeGen::make(const Eigen::Vector3d &root_pos, double trunk_radius, double random_factor)
{
  com.setZero();
  total_mass = 0.0;
  root_ = root_pos;

  Pose base(root_pos, Eigen::Quaterniond(Eigen::AngleAxisd(random_factor * random(0.0, 2.0 * kPi), Eigen::Vector3d(0, 0, 1))));
  Branch branch;
  branch.tip = root_pos;
  branch.parent_index = -1;
  branch.radius = trunk_radius;
  branches_.push_back(branch);
  addBranch(0, base, trunk_radius, random_factor);

  com /= total_mass;
  // std::cout << "COM: " << COM.transpose() << ", grad = " << COM[2]/trunkRadius << std::endl;
  double scale = branchGradient / (com[2] / trunk_radius);
  for (auto &leaf : leaves_) 
    leaf = root_pos + (leaf - root_pos) * scale;
  for (auto &start : ray_starts_) 
    start = root_pos + (start - root_pos) * scale;
  for (auto &end : ray_ends_) 
    end = root_pos + (end - root_pos) * scale;
  for (auto &branch : branches_) 
  {
    branch.tip = root_pos + (branch.tip - root_pos) * scale;
    branch.radius *= scale;
  }
}

bool TreeGen::makeFromString(const std::string &line)
{
  int num_commas = (int)std::count(line.begin(), line.end(), ',');
  if ((num_commas % 5) != 4) // badly formatted
    return false; 
  int num_sections = (num_commas + 1)/5;
  std::istringstream ss(line);
  for (int s = 0; s<num_sections; s++)
  {
    Branch section;
    std::string token;
    for (int i = 0; i<3; i++)
    {
      std::getline(ss, token, ',');
      section.tip[i] = std::stod(token.c_str());
    }
    if (s == 0)
    {
      root_ = section.tip - Eigen::Vector3d(0,0,0.1);
    }
    std::getline(ss, token, ',');
    section.radius = std::stod(token.c_str());
    std::getline(ss, token, ',');
    section.parent_index = std::stoi(token.c_str());
    if (branches_.size() > 0 && branches_.back().tip == section.tip)
    {
      std::cout << "Error: zero length branch at " << s << std::endl;
      return false;
    }
    branches_.push_back(section);
  }
  return true;
}

// create a set of rays covering the tree at a roughly uniform distribution
void TreeGen::generateRays(double ray_density)
{
  ASSERT(branches_.size() > 0);
  const double path_trunk_multiplier = 12.0; // observe the tree from this many trunk radii away
  const double ground_path_multiplier = 5.0; // observe the tree from this many trunk radii in height
  const double flight_path_multiplier = 20.0; // overhead path is at this many trunk radii above the ground
  double path_radius = branches_[0].radius * path_trunk_multiplier;
  double ring_heights[2] = {branches_[0].radius * ground_path_multiplier, branches_[0].radius * flight_path_multiplier};
  Eigen::Vector3d root = branches_[0].tip;

  std::vector<double> cumulative_size(branches_.size());
  cumulative_size[0] = 0;
  for (int i = 1; i < (int)branches_.size(); i++)
  {
    Branch &branch = branches_[i];
    Branch &parent_branch = branches_[branch.parent_index];
    double area = (branch.tip - parent_branch.tip).norm() * 2.0 * kPi * (branch.radius + parent_branch.radius) /
                  2.0;  // slightly approximate
    cumulative_size[i] = cumulative_size[i - 1] + area;
  }

  int num_rays = (int)(ray_density * cumulative_size.back());
  double area_per_ray = cumulative_size.back() / (double)num_rays;
  double total_area = 0.0;
  int i = 0;
  for (int j = 0; j < num_rays; j++)
  {
    total_area += area_per_ray;
    while (i < (int)cumulative_size.size() - 1 && cumulative_size[i] < total_area) i++;

    Branch &branch = branches_[i];
    Branch &parent_branch = branches_[branch.parent_index];

    // simplest is to randomise a point on a cone.... it will have overlap and gaps, but it is a starting point...
    double t = random(0.0, 1.0);
    double r = branch.radius + (parent_branch.radius - branch.radius) * t;
    Eigen::Vector3d online = branch.tip + (parent_branch.tip - branch.tip) * t;
    double angle = random(0.0, 2.0 * kPi);
    Eigen::Vector3d up = (branch.tip - parent_branch.tip).normalized();
    Eigen::Vector3d side = up.cross(Eigen::Vector3d(1, 2, 3));
    Eigen::Vector3d fwd = up.cross(side).normalized();
    side = fwd.cross(up);
    Eigen::Vector3d offset = side * sin(angle) + fwd * cos(angle);
    Eigen::Vector3d pos = online + offset * r;
    ray_ends_.push_back(pos);
    Eigen::Vector3d from = Eigen::Vector3d(random(-1, 1), random(-1, 1), random(-1, 1));
    if (from.dot(offset) < 0.0)
      from = -from;
    
    // closest point to branch position on circle, that is not passing through that branch. 
    // This circles are similar to if the tree was scanned by a scanner moving on two circular paths
    // Note that an ideal ray start would never pass through any other branch, but this is too expensive for a simple
    // example tree. 
    Eigen::Vector3d best_start(0,0,0);
    double min_dist2 = 0;
    for (int k = 0; k<2; k++)
    {
      Eigen::Vector3d root2 = root + Eigen::Vector3d(0,0,ring_heights[k]);
      Eigen::Vector3d to_path = pos - root;
      to_path[2] = 0.0;
      to_path.normalize();
      double radius = path_radius * (k==0 ? 0.4 : 1.0);
      Eigen::Vector3d start = root2 + to_path * radius;
      if ((start - pos).dot(offset) < 0.0)
      {
        Eigen::Vector3d side(to_path[1], -to_path[2], 0.0);
        if ((pos - root2).dot(side) > 0.0)
          start = root2 + side * radius;
        else 
          start = root2 - side * radius;
      }
      double dist2 = (start - pos).squaredNorm();
      if (k == 0 || dist2 < min_dist2)
      {
        best_start = start;
        min_dist2 = dist2;
      }
    }
    ray_starts_.push_back(best_start);
  }
}
} // ray