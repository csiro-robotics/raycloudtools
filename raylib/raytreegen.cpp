// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytreegen.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>

namespace ray
{
// given the input @c main_branch_angle of the larger branch, we can obtain the second branch angle and the
// two branch radii (for a base radius of 1). These branch attributes are unique if the branch follows
// Leonardo's rule (the cross sectional area is unchanged by the bifurcation) and the rule that the centre of
// mass direction is unchanged by the bifurcation.
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
/// returns the branch angle that best fits a given @c tip_angle
double getMainBranchAngle(double tip_angle)
{
  double x = tip_angle * (double)kBranchAngleSize / (0.5 * kPi);
  int i = std::max(0, std::min((int)x, kBranchAngleSize - 2));
  double blend = std::max(0.0, std::min(x - (double)i, 1.0));
  return branch_angle_lookup[i] * (1.0 - blend) + branch_angle_lookup[i + 1] * blend;
}

/// pre-fill the lookup table that approximates the choice of branch angles
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

void TreeGen::addBranch(int parent_index, Pose pose, double radius, const TreeParams &params)
{
  if (radius < params.min_branch_radius)
  {
    leaves_.push_back(pose.position);
    return;
  }
  Eigen::Vector3d p1 = pose.position;
  double rand_scale = random(1.0 - params.random_factor, 1.0 + params.random_factor);
  pose.position = pose * Eigen::Vector3d(0, 0, radius * branchGradient * rand_scale);
  double phi = (sqrt(5) + 1.0) / 2.0;
  rand_scale = random(1.0 - params.random_factor, 1.0 + params.random_factor);
  pose.rotation =
    pose.rotation * Eigen::Quaterniond(Eigen::AngleAxisd(2.0 * kPi * phi * rand_scale, Eigen::Vector3d(0, 0, 1)));
  Segment branch;
  branch.tip = pose.position;
  branch.radius = radius;
  branch.parent_id = parent_index;
  int index = int(segments_.size());
  segments_.push_back(branch);
  com += radius * radius * radius * (branch.tip + p1) / 2.0;
  total_mass += radius * radius * radius;

  Pose child1 = pose, child2 = pose;
  double angle1, angle2, radius1, radius2;
  rand_scale = random(1.0 - params.random_factor, 1.0 + params.random_factor);
  angle1 = getMainBranchAngle(kPi * 0.5 - pitchAngle * rand_scale);  // splitAngle*0.5 * mainBranchAngleRatio;
  getBranchInfo(angle1, angle2, radius1, radius2);
  rand_scale = random(1.0 - params.random_factor, 1.0 + params.random_factor);
  angle2 *= rand_scale;
  radius1 *= radius;
  radius2 *= radius;

  Eigen::Quaterniond q1(Eigen::AngleAxisd(angle1, Eigen::Vector3d(1, 0, 0)));
  Eigen::Quaterniond q2(Eigen::AngleAxisd(angle2, Eigen::Vector3d(-1, 0, 0)));

  child1.rotation = child1.rotation * q1;
  child2.rotation = child2.rotation * q2;
  addBranch(index, child1, radius1, params);
  addBranch(index, child2, radius2, params);
}

// create the tree structure, and list of leaf points
void TreeGen::make(const TreeParams &params)
{
  com.setZero();
  total_mass = 0.0;
  Pose base(segments_[0].tip, Eigen::Quaterniond(Eigen::AngleAxisd(params.random_factor * random(0.0, 2.0 * kPi),
                                                                   Eigen::Vector3d(0, 0, 1))));

  addBranch(0, base, segments_[0].radius, params);

  com /= total_mass;
  // having made the tree, we now scale the whole thing in order that it matches the expected tapering gradient
  double scale = branchGradient / (com[2] / segments_[0].radius);
  const Eigen::Vector3d &root = TreeStructure::root();
  for (auto &leaf : leaves_) leaf = root + (leaf - root) * scale;
  for (auto &start : ray_starts_) start = root + (start - root) * scale;
  for (auto &end : ray_ends_) end = root + (end - root) * scale;
  for (auto &branch : segments_)
  {
    branch.tip = root + (branch.tip - root) * scale;
    branch.radius *= scale;
  }
}

// create a set of rays covering the tree at a roughly uniform distribution
void TreeGen::generateRays(double ray_density)
{
  ASSERT(segments_.size() > 0);
  const double path_trunk_multiplier = 12.0;   // observe the tree from this many trunk radii away
  const double ground_path_multiplier = 5.0;   // observe the tree from this many trunk radii in height
  const double flight_path_multiplier = 20.0;  // overhead path is at this many trunk radii above the ground
  double path_radius = segments_[0].radius * path_trunk_multiplier;
  double ring_heights[2] = { segments_[0].radius * ground_path_multiplier,
                             segments_[0].radius * flight_path_multiplier };
  Eigen::Vector3d root = segments_[0].tip;

  std::vector<double> cumulative_size(segments_.size());
  cumulative_size[0] = 0;
  for (int i = 1; i < (int)segments_.size(); i++)
  {
    Segment &branch = segments_[i];
    Segment &parent_branch = segments_[branch.parent_id];
    double area = (branch.tip - parent_branch.tip).norm() * 2.0 * kPi * (branch.radius + parent_branch.radius) /
                  2.0;  // slightly approximate
    area *= random(0.10, 1.0);
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

    Segment &branch = segments_[i];
    Segment &parent_branch = segments_[branch.parent_id];

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
    {
      from = -from;
    }

    // closest point to branch position on circle, that is not passing through that branch.
    // This circles are similar to if the tree was scanned by a scanner moving on two circular paths
    // Note that an ideal ray start would never pass through any other branch, but this is too expensive for a simple
    // example tree.
    Eigen::Vector3d best_start(0, 0, 0);
    double min_dist2 = 0;
    for (int k = 0; k < 2; k++)
    {
      Eigen::Vector3d root2 = root + Eigen::Vector3d(0, 0, ring_heights[k]);
      Eigen::Vector3d to_path = pos - root;
      to_path[2] = 0.0;
      to_path.normalize();
      double radius = path_radius * (k == 0 ? 0.4 : 1.0);
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
}  // namespace ray