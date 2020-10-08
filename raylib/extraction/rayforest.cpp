// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

namespace ray
{

void Forest::searchTrees(const std::vector<TreeNode> &trees, int ind, double error, double length_per_radius, double ground_height, std::vector<int> &indices)
{
  if (trees[ind].children[0] == -1)
  {
    if (trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_))
      indices.push_back(ind);
    return;
  }
  int ind0 = trees[ind].children[0];
  double base0 = trees[ind0].height() - length_per_radius * trees[ind0].crownRadius();
  double error0 = abs(base0 - ground_height);
  int ind1 = trees[ind].children[1];
  double base1 = trees[ind1].height() - length_per_radius * trees[ind1].crownRadius();
  double error1 = abs(base1 - ground_height);
      
  if (error < std::min(error0, error1) && trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_)) // we've found the closest, so end loop
  {
    indices.push_back(ind);
    return;
  }
  searchTrees(trees, ind0, error0, length_per_radius, ground_height, indices);
  searchTrees(trees, ind1, error1, length_per_radius, ground_height, indices);
}

struct Point 
{ 
  int x, y, index; // if index == -2 then we are merging
  double height;
};
struct PointCmp 
{
  bool operator()(const Point& lhs, const Point& rhs) const 
  { 
    return lhs.height > rhs.height; 
  }
};

// extract the ray cloud canopy to a height field, then call the heightfield based forest extraction
void Forest::extract(const Cloud &cloud)
{
  cloud.calcBounds(&min_bounds_, &max_bounds_);
  double voxel_width = 4.0 * cloud.estimatePointSpacing();

  double width = (max_bounds_[0] - min_bounds_[0])/voxel_width;
  double length = (max_bounds_[1] - min_bounds_[1])/voxel_width;
  Eigen::Vector2i grid_dims(ceil(width), ceil(length));
  Eigen::ArrayXXd highs = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], -1e10);
  Eigen::ArrayXXd lows = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], 1e10);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    const Eigen::Vector3d &p = cloud.ends[i];
    Eigen::Vector3d pos = (p- min_bounds_)/voxel_width;
    double &h = highs((int)pos[0], (int)pos[1]);
    h = std::max(h, p[2]);
    double &l = lows((int)pos[0], (int)pos[1]);
    l = std::min(l, p[2]);
  }

  extract(highs, lows, voxel_width);
}

// Extraction uses a hierarchical watershed algorithm:
// 1. find all the highest points, give them unique labels, sort them from highest to lowest
// 2. for each highest point, make it black, then give all the non-black neighbours the same label and 
//    add them to the sorted list
void Forest::extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, double voxel_width)
{
  voxel_width_ = voxel_width;
  heightfield_ = highs; 
  lowfield_ = lows; 
  lowest_point_ = lows.minCoeff();

  // so we have some ground in lowfield_ that is beneath trees, we also have some of heightfield_ that is treeless ground.
  // The first we use directly in the ground estimation, the second we have to estimate based on its curvature and location...
  // If it is too planar, or the tree base is way too low compared to the others, then it is probably ground.

  indexfield_ = Eigen::ArrayXXi::Constant(heightfield_.rows(), heightfield_.cols(), -1);

  std::vector<TreeNode> trees;
  std::set<int> heads;
  hierarchicalWatershed(trees, heads);

  std::cout << "number of raw candidates: " << trees.size() << " number largest size: " << heads.size() << std::endl;

  calculateTreeParaboloids(trees);
  drawLowfield("lowfield.png", trees);
  drawSegmentation("segmented.png", trees);

  double ground_height = estimateRoundnessAndGroundHeight(trees);
  /***************************************************************************************/
  // The last step is to take the ground surface (currently flat: 'b') and scale estimation 'a'
  // and use this to estimate the set of tree heights
  // Note: we may still choose an extrapolation technique later... this is worth pursuing at some point
  std::vector<int> indices;
  for (auto &head: heads)
  {
    int ind = head;
    double base = trees[ind].height() - (1.0/tree_roundness) * trees[ind].crownRadius();
    double error = abs(base - ground_height);
    searchTrees(trees, head, error, 1.0/tree_roundness, ground_height, indices);
  }
  for (auto &ind: indices)
    result_.tree_tips.push_back(trees[ind].tip());
  result_.ground_height = ground_height;
  result_.treelength_per_crownradius = 1.0/tree_roundness;

  drawTrees("result_trees.png", result_);
}

double Forest::estimateRoundnessAndGroundHeight(std::vector<TreeNode> &trees)
{
  const double radius_to_height = 1.0;
  // Now collate weighted height vs crown radius data
  std::vector<Vector4d> data;
  // draw scatter plot
  double max_tree_height = -1e10;
  double min_tree_height = 1e10;
  double max_strength = 0.0;
  const bool sqrtWeight = false; // weights by the sqrt of the area... this gives more weigth to the small crowns
  for (auto &tree: trees)
  {
    double area = tree.numPoints() * ray::sqr(voxel_width_);
    double strength = sqrtWeight ? std::sqrt(area) : area;
    double ratio = tree.sum_square_residual/(1e-10 + tree.sum_square_total);
    double R2 = 1.0 - std::sqrt(ray::clamped(ratio, 0.0, 1.0)); 
    strength *= R2; // TODO: whether this helps is dubious, and it is weird that it gives values outside 0-1... more testing needed.
    max_strength = std::max(max_strength, strength);

    if (tree.validParaboloid(max_tree_canopy_width, voxel_width_))
    {
      max_tree_height = std::max(max_tree_height, tree.height());
      min_tree_height = std::min(min_tree_height, tree.height());
      data.push_back(Vector4d(tree.height(), tree.crownRadius(), strength, 1));
    }
  }
  // now add the under-canopy possible ground points:
  for (int i = 0; i<lowfield_.rows(); i++)
  {
    for (int j = 0; j<lowfield_.cols(); j++)
    {
      double weight = sqrtWeight ? voxel_width_ : voxel_width_ * voxel_width_;
      int ind = indexfield_(i,j);
      // Three indicators of ground:
      // 1. long way below highest point
      if (lowfield_(i,j) < heightfield_(i,j)-min_ground_to_canopy_distance)
        data.push_back(Vector4d(lowfield_(i,j), 0.0, weight, 2));
      // 2. unclassified (too big a drop from the tree edges)
      else if (ind == -1)
        data.push_back(Vector4d(lowfield_(i,j), 0.0, weight, 3));
      // 3. classified, but the paraboloid is too flat, or sloped or concave
      else
      {
        while (trees[ind].attaches_to != -1)
          ind = trees[ind].attaches_to;
        if (!trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_)) // invalid tree shape, so assume it is ground
          data.push_back(Vector4d(lowfield_(i,j), 0.0, weight, 4));
      }
    }
  }

  // now analyse the data to get line of best fit:
  // line = radius = a * treetop_pos + b
  double a = tree_roundness;
  double b = 0.0; 
  double min_height_resolving = 0.5; // unlikely to get better resolution than this
  for (int it = 0; it<13; it++) // from mean to median to mode (ish)
  {
    double power = (double)it / 9.0; // just 1 is median... its very similar, but not obviously better
    Vector4d mean(0,0,0,0);
    double total_weight = 0.0;
    for (auto &point: data)
    {
      double error = abs(point[1] - (a * point[0] + b));
      double weight = point[2] / (min_height_resolving + pow(error, power));
      mean += point * weight;
      total_weight += weight;
    }
    mean /= total_weight;
    if (tree_roundness == 0.0 && average_height == 0.0) // do line of best fit:
    {
      double total_xy = 0.0;
      double total_xx = 0.0;
      for (auto &point: data)
      {
        double error = abs(point[1] - (a * point[0] + b));
        double weight = point[2] / (min_height_resolving + pow(error, power));
        Vector4d p = point - mean;
        total_xy += p[0]*p[1] * weight;
        total_xx += p[0]*p[0] * weight;
      }
      a = total_xy / std::max(1e-10, total_xx);
      b = mean[1] - a*mean[0];

      // if the estimated ground is a lot higher than the lowest recorded point then limit its height
      double height = -b/a;
      const double ground_window = 3.0;
      if (height > lowest_point_ + ground_window)
        b = -a*(lowest_point_ + ground_window);
    }
    else if (average_height < 0.0) // in this case we force the ground height to be the lowest point
    {
      if (tree_roundness == 0.0)
        a = mean[1] / (mean[0] - lowest_point_);
      b = -lowest_point_ * a;
    }
    else 
    { 
      if (tree_roundness == 0.0)
        a = mean[1] / average_height;
      else if (average_height == 0.0)
        average_height = mean[1]/a;
      b = (average_height - mean[0]) * a;
    }
  }
  drawGraph("graph_curv.png", data, min_tree_height, max_tree_height, 20.0 * ray::sqr(radius_to_height), max_strength, a, b);
  std::cout << "a: " << a << ", b: " << b << std::endl;
  a = 0.32; // actual value
  tree_roundness = a;
  double ground_height = -b/a;
  return ground_height;
}

void Forest::hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads)
{
  std::set<Point, PointCmp> basins;
  // 1. find highest points
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
      // Moore neighbourhood
      double height = heightfield_(x, y);
      double max_h = 0.0;
      for (int i = std::max(0, x-1); i<= std::min(x+1, (int)heightfield_.rows()-1); i++)
        for (int j = std::max(0, y-1); j<= std::min(y+1, (int)heightfield_.cols()-1); j++)
          if (!(i==x && j==y))
            max_h = std::max(max_h, heightfield_(i, j));
      if (height > max_h)
      {
        Point p;
        p.x = x; p.y = y; p.height = height;
        p.index = (int)basins.size();
        basins.insert(p);
        heads.insert((int)trees.size());
        trees.push_back(TreeNode(x, y, height, voxel_width_));
      }
    }
  }
  std::cout << "initial number of peaks: " << trees.size() << std::endl;
  // now iterate until basins is empty
  int cnt = 0;
  // Below, don't divide by voxel_width, if you want to verify voxel_width independence
  int max_tree_pixel_width = (int)(max_tree_canopy_width / (double)voxel_width_); 
  while (!basins.empty())
  {
    Point p = *basins.begin();
    int x = p.x;
    int y = p.y;
    basins.erase(p); // removes it from basins. p still exists

    if (p.index == -2) // a merge request
    {
      int p_head = x;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
      int q_head = y;
      while (trees[q_head].attaches_to != -1)
        q_head = trees[q_head].attaches_to;
      if (p_head != q_head)
      {
        TreeNode &p_tree = trees[p_head];
        TreeNode &q_tree = trees[q_head];
        Eigen::Vector2i mx = ray::maxVector2(p_tree.max_bound, q_tree.max_bound);
        Eigen::Vector2i mn = ray::minVector2(p_tree.min_bound, q_tree.min_bound);
        mx -= mn;
        if (std::max(mx[0], mx[1]) <= max_tree_pixel_width)
        {
          int new_index = (int)trees.size();
          TreeNode node;
          node.curv_mat = p_tree.curv_mat + q_tree.curv_mat;
          node.curv_vec = p_tree.curv_vec + q_tree.curv_vec;
          node.min_bound = p_tree.min_bound;
          node.max_bound = p_tree.max_bound;
          node.updateBound(q_tree.min_bound, q_tree.max_bound);
          node.children[0] = p_head;  
          node.children[1] = q_head;
          node.abcd = node.curv_mat.ldlt().solve(node.curv_vec);

          if (node.validParaboloid(max_tree_canopy_width, voxel_width_)) 
          {
            heads.erase(p_head);
            heads.erase(q_head);
            heads.insert(new_index);
            p_tree.attaches_to = new_index;
            q_tree.attaches_to = new_index;
            trees.push_back(node); // danger, this can invalidate the p_tree reference
          }
        }
      }
      continue;
    }    
    indexfield_(x, y) = p.index;

    int xs[4] = {x-1, x, x, x+1};
    int ys[4] = {y, y+1, y-1, y};
    for (int i = 0; i<4; i++)
    {
      if (xs[i] < 0 || xs[i] >= indexfield_.rows())
        continue;
      if (ys[i] < 0 || ys[i] >= indexfield_.cols())
        continue;
      int p_head = p.index;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
        
      int xx = xs[i];
      int yy = ys[i];
      int &ind = indexfield_(xx, yy);

      int q_head = ind;
      if (q_head != -1)
      {
        while (trees[q_head].attaches_to != -1)
          q_head = trees[q_head].attaches_to;
      }

      if (ind != -1 && p_head != q_head)
      {
        TreeNode &p_tree = trees[p_head];
        TreeNode &q_tree = trees[q_head];
        cnt++;
//        if (verbose && !(cnt%50)) // I need a way to visualise the hierarchy here!
//          drawSegmentation("segmenting.png", trees);
        Eigen::Vector2i mx = ray::maxVector2(p_tree.max_bound, q_tree.max_bound);
        Eigen::Vector2i mn = ray::minVector2(p_tree.min_bound, q_tree.min_bound);
        mx -= mn;
        bool merge = std::max(mx[0], mx[1]) <= max_tree_pixel_width;
        if (merge)
        {
          const double flood_merge_scale = 2.0; // 1 merges immediately, infinity never merges
          // add a merge task:
          Eigen::Vector2d mid = Eigen::Vector2d(xx, yy) * voxel_width_;
          double p_sqr = (Eigen::Vector2d(p_tree.peak[0], p_tree.peak[1]) - mid).squaredNorm();
          double q_sqr = (Eigen::Vector2d(q_tree.peak[0], q_tree.peak[1]) - mid).squaredNorm();
          double blend = p_sqr / (p_sqr + q_sqr);
          double flood_base = p_tree.peak[2]*(1.0-blend) + q_tree.peak[2]*blend;
          double low_flood_height = flood_base - p.height;

          Point q;
          q.x = p_head; q.y = q_head; 
          q.index = -2;
          q.height = flood_base - low_flood_height * flood_merge_scale;
          basins.insert(q);
        }
      }
      if (ind == -1 && heightfield_(xx, yy) > -1e10) 
      {
        Point q;
        q.x = xx; q.y = yy; q.index = p.index;
        q.height = heightfield_(xx, yy);
        if ((p.height - q.height) < 4.0)
        {
          ind = p.index;
          basins.insert(q);
          trees[p_head].addSample(xx*voxel_width_, yy*voxel_width_, q.height);
          trees[p_head].updateBound(Eigen::Vector2i(xx, yy), Eigen::Vector2i(xx, yy));
        }
      }
    }
  }
}

void Forest::calculateTreeParaboloids(std::vector<TreeNode> &trees)
{
  // first I'll have to solve each tree parabola and store it in a vector:
  for (auto &tree: trees)
    tree.abcd = tree.curv_mat.ldlt().solve(tree.curv_vec);

  // for each pixel, we have to update accuracy data on each tree.
  for (int i = 0; i < indexfield_.rows(); i++)
  {
    for (int j = 0; j < indexfield_.cols(); j++)
    {
      int ind = indexfield_(i, j);
      double x = i*voxel_width_;
      double y = j*voxel_width_;
      if (ind != -1)
      {
        trees[ind].sum_square_residual += ray::sqr(heightfield_(i,j) - trees[ind].heightAt(x,y));
        trees[ind].sum_square_total += ray::sqr(heightfield_(i,j) - trees[ind].avgHeight());
        while (trees[ind].attaches_to != -1)
        {
          ind = trees[ind].attaches_to;
          trees[ind].sum_square_residual += ray::sqr(heightfield_(i,j) - trees[ind].heightAt(x,y));
          trees[ind].sum_square_total += ray::sqr(heightfield_(i,j) - trees[ind].avgHeight());
        }
      }
    }
  }
}
}