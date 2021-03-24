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
#include <queue>

namespace ray
{

double Forest::searchTrees(const std::vector<TreeNode> &trees, int ind, double length_per_radius, std::vector<int> &indices)
{
  double base = trees[ind].node.height() - length_per_radius * trees[ind].node.crownRadius();
  double error = abs(base - trees[ind].ground_height);
  if (trees[ind].children[0] == -1)
  {
    if (trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_))
    {
      indices.push_back(ind);
      return error;
    }
    return 1e20;
  }
  std::vector<int> child_indices[2];
  int ind0 = trees[ind].children[0];
  int ind1 = trees[ind].children[1];
  double child_error = searchTrees2(trees, ind0, length_per_radius, child_indices[0]);
  if (ind1 != -1)
  {
    child_error = (child_error + searchTrees2(trees, ind1, length_per_radius, child_indices[1])) / 2.0; // mean error
  }

  if (error < child_error && trees[ind].validParaboloid(max_tree_canopy_width, voxel_width_))
  {
    indices.push_back(ind);
    return error;
  }

  indices.insert(indices.end(), child_indices[0].begin(), child_indices[0].end());
  indices.insert(indices.end(), child_indices[1].begin(), child_indices[1].end());
  return child_error;
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
    return lhs.height < rhs.height; 
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
      if (heightfield_(x, y) < lowfield_(x, y)+undercroft_height)
      {
        heightfield_(x, y) = -1e10;
        count++;
      }
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
    searchTrees(trees, head, 1.0/tree_roundness, indices);
  }
  for (auto &ind: indices)
  {
    Result result;
    result.tree_tip = trees[ind].node.mean(); // peak;//tip();
    result.tree_tip[2] = trees[ind].peak[2];
    int x = (int) result.tree_tip[0];
    int y = (int) result.tree_tip[1];
    if (x < 0 || x >= lowfield_.rows())
      continue;
    if (y < 0 || y >= lowfield_.cols())
      continue;
    result.ground_height = lowfield_(x, y);
    result.radius = trees[ind].approx_radius;
    result.curvature = trees[ind].node.curvature();
    results_.push_back(result);
  }

  drawTrees("result_trees.png", results_, (int)heightfield_.rows(), (int)heightfield_.cols());

  std::vector<TreeNode> selected;
  for (auto &ind: indices)
    selected.push_back(trees[ind]);
    
  drawFinalSegmentation("result_tree_shapes.png", trees, indices);
}

void Forest::hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads)
{
  std::priority_queue<Point, std::vector<Point>, PointCmp> basins;
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
        basins.push(p);
        heads.insert((int)trees.size());
        indexfield_(x, y) = p.index;
        if (trees.size() == 0)
          std::cout << "x: " << x << ",  y: " << y << ", index field: " << indexfield_(x, y) << std::endl;      
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
    Point p = basins.top();
    basins.pop(); // removes it from basins. p still exists
    int x = p.x;
    int y = p.y;

    if (p.index == -2) // a merge request
    {
      int p_head = x;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
      int q_head = y;
      while (trees[q_head].attaches_to != -1)
        q_head = trees[q_head].attaches_to;
      if (p_head != q_head) // already merged
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
          node.min_bound = p_tree.min_bound;
          node.max_bound = p_tree.max_bound;
          node.updateBound(q_tree.min_bound, q_tree.max_bound);
          node.children[0] = p_head;  
          node.children[1] = q_head;
          node.peak = p_tree.peak[2] > q_tree.peak[2] ? p_tree.peak : q_tree.peak;

 //         if (node.validParaboloid(max_tree_canopy_width, voxel_width_)) 
          {
            heads.erase(p_head);
            heads.erase(q_head);
            heads.insert(new_index);
            p_tree.attaches_to = new_index;
            q_tree.attaches_to = new_index;
            trees.push_back(node); // danger, this can invalidate the p_tree reference

            // Below: no noticeable difference in quality, so leaving out
            #if 0 // adds in a node just at the merge point, for more fidelity
            trees.back().attaches_to = (int)trees.size();
            TreeNode node2 = node;
            node.children[0] = trees.size()-1;
            node.children[1] = -1;
            heads.erase(new_index);
            heads.insert(trees.size());
            trees.push_back(node2);
            #endif
          }

        }
      }
      continue;
    }    

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
        if (merge)
        {
          const double flood_merge_scale = 2.0; // 1 merges immediately, infinity never merges
          // add a merge task:
          Eigen::Vector2d mid = Eigen::Vector2d(xx, yy) * voxel_width_;
          Eigen::Vector2d ptree(p_tree.peak[0], p_tree.peak[1]);
          Eigen::Vector2d qtree(q_tree.peak[0], q_tree.peak[1]);
          double blend = (mid - ptree).dot(qtree-ptree) / (qtree-ptree).squaredNorm();
          double flood_base = p_tree.peak[2]*(1.0-blend) + q_tree.peak[2]*blend;
          double low_flood_height = flood_base - p.height;

          Point q;
          q.x = p_head; q.y = q_head; 
          q.index = -2;
          q.height = flood_base - low_flood_height * flood_merge_scale;
          basins.push(q);
        }
      }
      if (ind == -1 && heightfield_(xx, yy) > -1e10) 
      {
        Point q;
        q.x = xx; q.y = yy; q.index = p.index;
        q.height = heightfield_(xx, yy);
        if ((p.height - q.height) < maximum_drop_within_tree)
        {
  /*        if (verbose && !(cnt%500)) // I need a way to visualise the hierarchy here!
          {
            drawSegmentation("segmenting.png", trees);
          }*/
          cnt++;
          ind = p.index;
          basins.push(q);
          trees[p_head].updateBound(Eigen::Vector2i(xx, yy), Eigen::Vector2i(xx, yy));        
        } 
      }
    }
  }
}

void Forest::calculateTreeParaboloids(std::vector<TreeNode> &trees)
{
  std::vector<std::vector<Eigen::Vector3d> > point_lists(trees.size());
  for (int x = 0; x<indexfield_.rows(); x++)
  {
    for (int y = 0; y<indexfield_.cols(); y++)
    {
      int ind = indexfield_(x, y);
      if (ind < 0)
        continue;
      while (ind >= 0)
      {
        point_lists[ind].push_back(Eigen::Vector3d((double)x + 0.5, (double)y + 0.5, heightfield_(x, y)));
        ind = trees[ind].attaches_to;
      }
    }
  }
  for (size_t i = 0; i<trees.size(); i++)
  {
    auto &tree = trees[i];
    tree.approx_radius = std::sqrt((double)point_lists[i].size() / kPi);
    int x = (int)(tree.peak[0]/voxel_width_);
    int y = (int)(tree.peak[1]/voxel_width_);   
    x = std::max(0, std::min(x, (int)lowfield_.rows()-1));
    y = std::max(0, std::min(y, (int)lowfield_.cols()-1));
    tree.ground_height = lowfield_(x, y); 
    TreeNode::Node node;
    for (auto &pt: point_lists[i])
      node.add(pt[0], pt[1], pt[2], 1);
    const int num_iterations = 10;
    for (int it = 1; it<num_iterations; it++)
    {
      node.abcd = node.curv_mat.ldlt().solve(node.curv_vec);
      node.curv_mat.setZero(); 
      node.curv_vec.setZero();
      for (auto &pt: point_lists[i])
      {
        double h = node.heightAt(pt[0], pt[1]);
        double error = h - pt[2];
        const double eps = 1e-2;
        node.add(pt[0], pt[1], pt[2], 1.0/std::max(eps, std::abs(error))); // 1/e reweighting gives a median paraboloid
      }
    }
    node.abcd = node.curv_mat.ldlt().solve(node.curv_vec);
    tree.node = node;
  }
}
}
