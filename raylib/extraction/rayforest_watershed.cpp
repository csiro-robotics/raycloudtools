// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>
#include "../rayply.h"

namespace ray
{
void Forest::renderWatershed(const std::string &cloud_name_stub, std::vector<TreeNode> &trees, std::set<int> &indices)
{
  if (!verbose)
    return;
  std::vector<Eigen::Vector3d> cloud_points;
  std::vector<double> times;
  std::vector<RGBA> colours;
  RGBA colour;
  colour.alpha = 255;

  for (int x = 0; x < indexfield_.rows(); x++)
  {
    for (int y = 0; y < indexfield_.cols(); y++)
    {
      int ind = indexfield_(x, y);
      if (ind == -1)
      {
        continue;
      }
      while (trees[ind].attaches_to != -1) ind = trees[ind].attaches_to;
      if (trees[ind].area < min_area_)
      {
        continue;
      }
      srand(1 + ind);
      colour.red = (uint8_t)(rand() % 256);
      colour.green = (uint8_t)(rand() % 256);
      colour.blue = (uint8_t)(rand() % 256);
      Eigen::Vector3d pos =
        min_bounds_ + voxel_width_ * Eigen::Vector3d(0.5 + static_cast<double>(x), 0.5 + static_cast<double>(y), 0);
      pos[2] = original_heightfield_(x, y);
      cloud_points.push_back(pos);
      times.push_back(0.0);
      colours.push_back(colour);
    }
  }

  for (auto &ind : indices)
  {
    Eigen::Vector3d tip;
    srand(1 + ind);
    colour.red = (uint8_t)(rand() % 256);
    colour.green = (uint8_t)(rand() % 256);
    colour.blue = (uint8_t)(rand() % 256);
    if (trees[ind].area < min_area_)
    {
      continue;
    }
    double z_max = 2.0;
    if (trees[ind].trunk_id >= 0)
    {
      z_max = 4.0;
    }
    if (findSpace(trees[ind], tip))
    {
      Eigen::Vector3d base = min_bounds_ + tip;
      base[2] = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
      const double rad = tip[2] / approx_height_per_radius_;

      for (double z = 0.0; z < z_max; z += 0.3)
      {
        for (double ang = 0.0; ang < 2.0 * kPi; ang += 0.3)
        {
          cloud_points.push_back(base + Eigen::Vector3d(rad * std::sin(ang), rad * std::cos(ang), z));
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }
  }

  // now add the space field:
  for (int i = 0; i < spacefield_.rows(); i++)
  {
    for (int j = 0; j < spacefield_.cols(); j++)
    {
      if (spacefield_(i, j) < 1.0)
      {
        const double height = lowfield_(i, j) + 0.2;
        const double x = min_bounds_[0] + static_cast<double>(i) * voxel_width_;
        const double y = min_bounds_[1] + static_cast<double>(j) * voxel_width_;
        cloud_points.push_back(Eigen::Vector3d(x, y, height));
        times.push_back(0.0);
        colour.red = colour.green = colour.blue = (uint8_t)(255.0 * spacefield_(i, j));
        colours.push_back(colour);
      }
    }
  }

  writePlyPointCloud(cloud_name_stub + "_watershed.ply", cloud_points, times, colours);
}

struct Point
{
  int x, y, index;  // if index == -2 then we are merging
  double height;
};
struct PointCmp
{
  bool operator()(const Point &lhs, const Point &rhs) const { return lhs.height < rhs.height; }
};

void Forest::hierarchicalWatershed(std::vector<TreeNode> &trees, std::set<int> &heads)
{
  // fast array lookup of trunk centres:
  Eigen::ArrayXXi trunkfield = Eigen::ArrayXXi::Constant(indexfield_.rows(), indexfield_.cols(), -1);
  for (int c = 0; c < static_cast<int>(trunks_.size()); c++)  // if there are known trunks, then include them...
  {
    auto &trunk = trunks_[c];
    const Eigen::Vector3i pos = ((trunk.first - min_bounds_) / voxel_width_).cast<int>();
    if (pos[0] < 0 || pos[0] >= trunkfield.rows() || pos[1] < 0 || pos[1] >= trunkfield.cols())
    {
      std::cout << "warning: trunk " << c << " location is out of bounds" << std::endl;
      continue;
    }
    trunkfield(pos[0], pos[1]) = c;
  }

  std::priority_queue<Point, std::vector<Point>, PointCmp> basins;
  // 1. find highest points
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
      // Moore neighbourhood
      const double height = heightfield_(x, y);
      double max_h = 0.0;
      for (int i = std::max(0, x - 1); i <= std::min(x + 1, static_cast<int>(heightfield_.rows()) - 1); i++)
      {
        for (int j = std::max(0, y - 1); j <= std::min(y + 1, static_cast<int>(heightfield_.cols()) - 1); j++)
        {
          if (!(i == x && j == y))
          {
            max_h = std::max(max_h, heightfield_(i, j));
          }
        }
      }
      if (height > max_h && height > std::numeric_limits<double>::lowest())
      {
        Point p;
        p.x = x;
        p.y = y;
        p.height = height;
        p.index = static_cast<int>(basins.size());
        basins.push(p);
        heads.insert(p.index);
        indexfield_(x, y) = p.index;
        trees.push_back(TreeNode(x, y, height, voxel_width_, trunkfield(x, y)));
      }
    }
  }

  std::cout << "initial number of peaks: " << trees.size() << std::endl;
  int cnt = 0;
  // now iterate until basins is empty
  // Below, don't divide by voxel_width, if you want to verify voxel_width independence
  int max_tree_pixel_width = static_cast<int>(max_tree_canopy_width / static_cast<double>(voxel_width_));
  while (!basins.empty())
  {
    const Point p = basins.top();
    basins.pop();  // removes it from basins. p still exists
    const int x = p.x;
    const int y = p.y;
    const int xs[4] = { x - 1, x, x, x + 1 };
    const int ys[4] = { y, y + 1, y - 1, y };
    for (int i = 0; i < 4; i++)
    {
      if (xs[i] < 0 || xs[i] >= indexfield_.rows())
      {
        continue;
      }
      if (ys[i] < 0 || ys[i] >= indexfield_.cols())
      {
        continue;
      }
      int p_head = p.index;
      while (trees[p_head].attaches_to != -1)
      {
        p_head = trees[p_head].attaches_to;
      }

      const int xx = xs[i];
      const int yy = ys[i];
      int &ind = indexfield_(xx, yy);

      int q_head = ind;
      if (q_head != -1)
      {
        while (trees[q_head].attaches_to != -1)
        {
          q_head = trees[q_head].attaches_to;
        }
      }

      if (ind != -1 && p_head != q_head)  // connecting separate trees, so trigger a future merge event
      {
        TreeNode &p_tree = trees[p_head];
        TreeNode &q_tree = trees[q_head];
        Eigen::Vector2i mx = p_tree.max_bound.cwiseMax(q_tree.max_bound);
        Eigen::Vector2i mn = p_tree.min_bound.cwiseMin(q_tree.min_bound);
        mx -= mn;

        const bool mergable = std::max(mx[0], mx[1]) <= max_tree_pixel_width;
        const bool separate_trunks = p_tree.trunk_id >= 0 && q_tree.trunk_id >= 0 && p_tree.trunk_id != q_tree.trunk_id;
        if (mergable && !separate_trunks)
        {
          // add a merge task:
          const Eigen::Vector2d mid = Eigen::Vector2d(xx, yy) * voxel_width_;
          const Eigen::Vector2d ptree(p_tree.peak[0], p_tree.peak[1]);
          const Eigen::Vector2d qtree(q_tree.peak[0], q_tree.peak[1]);
          const double blend =
            std::max(0.0, std::min((mid - ptree).dot(qtree - ptree) / (qtree - ptree).squaredNorm(), 1.0));
          const double peak = p_tree.peak[2] * (1.0 - blend) + q_tree.peak[2] * blend;
          const double drop = peak - p.height;
          const double tree_height = std::max(0.0, peak - lowfield_(xx, yy));

          // if there is no space under one of the two areas, then do the merge
          Eigen::Vector3d tip;
          const bool space_each = findSpace(p_tree, tip) && findSpace(q_tree, tip);

          const bool too_small = std::max(mx[0], mx[1]) <= 10;
          if (drop < tree_height * drop_ratio_ || too_small || !space_each)  // good to merge
          {
            const int new_index = static_cast<int>(trees.size());
            TreeNode node;
            node.peak = p_tree.peak[2] > q_tree.peak[2] ? p_tree.peak : q_tree.peak;
            node.min_bound = p_tree.min_bound;
            node.max_bound = p_tree.max_bound;
            node.updateBound(q_tree.min_bound, q_tree.max_bound);
            node.children[0] = p_head;
            node.children[1] = q_head;
            node.trunk_id = p_tree.trunk_id >= 0 ? p_tree.trunk_id : q_tree.trunk_id;
            node.height = p.height;

            heads.erase(p_head);
            heads.erase(q_head);
            heads.insert(new_index);
            p_tree.attaches_to = new_index;
            q_tree.attaches_to = new_index;
            trees.push_back(node);  // danger, this can invalidate the p_tree reference
          }
        }
      }
      if (ind == -1 && heightfield_(xx, yy) > std::numeric_limits<double>::lowest())  // adding a single pixel to a tree
      {
        Point q;
        q.x = xx;
        q.y = yy;
        q.height = heightfield_(xx, yy);
        cnt++;

        const int trunkid = trunkfield(xx, yy);
        if (trunkid >= 0)
        {
          if (trees[p_head].trunk_id == -1)
          {
            trees[p_head].trunk_id = trunkid;
          }
          else  // a second trunk on a downward slope, we'll have to make a whole new treenodde
          {
            p_head = static_cast<int>(trees.size());  // this new pixel will point to the new tree node here
            trees.push_back(TreeNode(xx, yy, q.height, voxel_width_, trunkid));
          }
        }
        q.index = p_head;
        ind = p_head;
        basins.push(q);
        trees[p_head].updateBound(Eigen::Vector2i(xx, yy), Eigen::Vector2i(xx, yy));
      }
    }
  }
}

}  // namespace ray
