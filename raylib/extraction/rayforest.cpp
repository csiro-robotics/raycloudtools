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

void Forest::searchTrees(const std::vector<TreeNode> &trees, int ind, double error, double length_per_radius, std::vector<int> &indices)
{
  if (trees[ind].children[0] == -1)
  {
    if (trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_))
      indices.push_back(ind);
    return;
  }
  int ind0 = trees[ind].children[0];
  double base0 = trees[ind0].height() - length_per_radius * trees[ind0].crownRadius();
  double error0 = abs(base0 - trees[ind0].ground_height);
  int ind1 = trees[ind].children[1];
  double base1 = trees[ind1].height() - length_per_radius * trees[ind1].crownRadius();
  double error1 = abs(base1 - trees[ind1].ground_height);
      
  if (error < std::min(error0, error1) && trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_)) // we've found the closest, so end loop
  {
    indices.push_back(ind);
    return;
  }
  searchTrees(trees, ind0, error0, length_per_radius, indices);
  searchTrees(trees, ind1, error1, length_per_radius, indices);
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
void Forest::extract(const Cloud &cloud, Mesh &mesh)
{
  cloud.calcBounds(&min_bounds_, &max_bounds_);
  double voxel_width = 1.0; // 4.0 * cloud.estimatePointSpacing();

  double width = (max_bounds_[0] - min_bounds_[0])/voxel_width;
  double length = (max_bounds_[1] - min_bounds_[1])/voxel_width;
  Eigen::Vector2i grid_dims(ceil(width), ceil(length));
  std::cout << "dimes for heigh: " << grid_dims.transpose() << std::endl;
  Eigen::ArrayXXd highs = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], -1e10);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    const Eigen::Vector3d &p = cloud.ends[i];
    Eigen::Vector3d pos = (p- min_bounds_)/voxel_width;
    double &h = highs((int)pos[0], (int)pos[1]);
    h = std::max(h, p[2]);
  }
  Eigen::ArrayXXd lows;
  mesh.toHeightField(lows, min_bounds_, max_bounds_, voxel_width);
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
  indexfield_ = Eigen::ArrayXXi::Constant(heightfield_.rows(), heightfield_.cols(), -1);
  // ignore the undercroft
  int count = 0;
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
 //     if (heightfield_(x, y) > -1e10)
 //       std::cout << "low: " << lowfield_(x, y) << ", high: " << heightfield_(x, y) << std::endl;
      if (heightfield_(x, y) < lowfield_(x, y)+undercroft_height)
      {
        heightfield_(x, y) = -1e10;
        count++;
      }
      int k = 200; // any larger and the bug appears
      if (x < 208-k || x > 224+k || y < 64-k || y > 81+k)
        heightfield_(x, y) = -1e10;
    }
  }
  std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;
  drawHeightField("highfield.png", heightfield_);
  drawHeightField("lowfield.png", lowfield_);
  std::vector<TreeNode> trees;
  std::set<int> heads;
  hierarchicalWatershed(trees, heads);

  std::cout << "number of raw candidates: " << trees.size() << " number largest size: " << heads.size() << std::endl;

  calculateTreeParaboloids(trees);
  drawSegmentation("segmented.png", trees);

  std::vector<int> indices;
  for (auto &head: heads)
  {
    int ind = head;
    double base = trees[ind].height() - (1.0/tree_roundness) * trees[ind].crownRadius();
    double error = abs(base - trees[ind].ground_height);
    searchTrees(trees, head, error, 1.0/tree_roundness, indices);
  }
  for (auto &ind: indices)
  {
    Result result;
    result.tree_tip = trees[ind].tip();
    result.ground_height = trees[ind].ground_height;
    results_.push_back(result);
  }

  drawTrees("result_trees.png", results_);
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
      if (height > max_h && height > -1e10)
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
  int cnt = 0;
  // now iterate until basins is empty
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

   //       if (node.validParaboloid(max_tree_canopy_width, voxel_width_)) 
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
        Eigen::Vector2i mx = ray::maxVector2(p_tree.max_bound, q_tree.max_bound);
        Eigen::Vector2i mn = ray::minVector2(p_tree.min_bound, q_tree.min_bound);
        mx -= mn;
        bool merge = std::max(mx[0], mx[1]) <= max_tree_pixel_width;
//        double mz = std::max(p_tree.peak[2], q_tree.peak[2]);
//        bool merge = std::max(mx[0], mx[1]) <= max_tree_pixel_width && (mz - p.height)<maximum_drop_within_tree;
        if (merge)
        {
          const double flood_merge_scale = 2.0; // 1 merges immediately, infinity never merges
          // add a merge task:
          Eigen::Vector2d mid = Eigen::Vector2d(xx, yy) * voxel_width_;
          double p_sqr = (Eigen::Vector2d(p_tree.peak[0], p_tree.peak[1]) - mid).squaredNorm();
          double q_sqr = (Eigen::Vector2d(q_tree.peak[0], q_tree.peak[1]) - mid).squaredNorm();
          double blend = p_sqr / (p_sqr + q_sqr);
          double flood_base = p_tree.peak[2]*(1.0-blend) + q_tree.peak[2]*blend;
//          double flood_base = std::max(p_tree.peak[2], q_tree.peak[2]); // p_tree.peak[2]*(1.0-blend) + q_tree.peak[2]*blend;
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
        if ((p.height - q.height) < maximum_drop_within_tree)
//        if (std::abs(p.height - q.height) < maximum_drop_within_tree)
        {
          if (verbose && !(cnt%500)) // I need a way to visualise the hierarchy here!
          {
            drawSegmentation("segmenting.png", trees);
            std::cout << "done" << std::endl;
          }
          cnt++;
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
  {
    tree.abcd = tree.curv_mat.ldlt().solve(tree.curv_vec);
    Eigen::Vector3d tip = tree.tip();
    int x = (int)((tip[0] - min_bounds_[0])/voxel_width_);
    int y = (int)((tip[1] - min_bounds_[1])/voxel_width_);   
    x = std::max(0, std::min(x, (int)lowfield_.rows()-1));
    y = std::max(0, std::min(y, (int)lowfield_.cols()-1));
    tree.ground_height = lowfield_(x, y); 
  }
}
}
