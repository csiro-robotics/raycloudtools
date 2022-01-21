// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"
#include "rayterrain.h"
#include "../rayconvexhull.h"
#include "../raymesh.h"
#include "../raycuboid.h"
#include "../rayply.h"
#include "../rayforestgen.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

namespace ray
{
static const double wrap_gradient = 1.0;

bool Forest::findSpace(const Cluster &cluster, const std::vector<Eigen::Vector3d> &points, Eigen::Vector3d &tip)
{
  if (cluster.ids.empty())
  {
    tip = (cluster.min_bound + cluster.max_bound) / 2.0;
    tip[2] = cluster.max_bound[2];
    return true; // TODO: do we believe trunks absolutely? or should we double check against the free space?
  }
  if (cluster.trunk_id >= 0) // if this cluster is associated with a trunk, then use the trunk location, not the centroid
  {
    tip = trunks_[cluster.trunk_id].first - min_bounds_;
    tip[2] = cluster.max_bound[2];    
    return true;
  }
  Eigen::Vector3d weighted_sum(0,0,0);
  double weight = 0.0;
  for (auto &i: cluster.ids)
  {
    weighted_sum += points[i][2] * points[i];
    weight += points[i][2];
  }
  tip = weighted_sum/weight;
  Eigen::Vector3d tip_local = tip/voxel_width_;
  const double search_down_gradient = 0.2;
  double radius = tip_local[2]*search_down_gradient;

  // now find the closest bit of space to put the tree in:
  int min_x = std::max(0, (int)(tip_local[0] - radius));
  int max_x = std::min((int)spacefield_.rows()-1, (int)(tip_local[0] + radius));
  int min_y = std::max(0, (int)(tip_local[1] - radius));
  int max_y = std::min((int)spacefield_.cols()-1, (int)(tip_local[1] + radius));
  double best_score = -1e10;
  int best_x = -1;
  int best_y = -1;
  for (int x = min_x; x<=max_x; x++)
  {
    for (int y = min_y; y<=max_y; y++)
    {
      double dist2 = sqr(((double)x-tip_local[0])/radius) + sqr(((double)y-tip_local[1])/radius);
      double score = spacefield_(x, y) - 0.25*dist2; // slight preference for result near the centroid
      if (score > best_score)
      {
        best_score = score;
        best_x = x;
        best_y = y;
      }
    }
  }
  if (best_score > 0.0) 
  {
    tip[0] = ((double)best_x+0.5)*voxel_width_;
    tip[1] = ((double)best_y+0.5)*voxel_width_;
    return true;
  }
  return false;
}

bool Forest::findSpace2(const TreeNode &node, Eigen::Vector3d &tip)
{
  bool calculate = true;
  Eigen::Vector2d mid = voxel_width_ * (node.min_bound + node.max_bound).cast<double>()/2.0;
  tip[0] = mid[0];
  tip[1] = mid[1];

  if (node.trunk_id >= 0) // if this node is associated with a trunk, then use the trunk location, not the centroid
  {
    tip = trunks_[node.trunk_id].first - min_bounds_;
    calculate = false;
  }

  Eigen::Vector3d tip_local = tip/voxel_width_;
  tip[2] = node.peak[2] - lowfield_((int)tip_local[0], (int)tip_local[1]);
  if (!calculate)
  {
    return true;
  }
  const double search_down_gradient = 0.2;
  double radius = tip_local[2]*search_down_gradient;

  // now find the closest bit of space to put the tree in:
  int min_x = std::max(0, (int)(tip_local[0] - radius));
  int max_x = std::min((int)spacefield_.rows()-1, (int)(tip_local[0] + radius));
  int min_y = std::max(0, (int)(tip_local[1] - radius));
  int max_y = std::min((int)spacefield_.cols()-1, (int)(tip_local[1] + radius));
  double best_score = -1e10;
  int best_x = -1;
  int best_y = -1;
  for (int x = min_x; x<=max_x; x++)
  {
    for (int y = min_y; y<=max_y; y++)
    {
      double dist2 = sqr(((double)x-tip_local[0])/radius) + sqr(((double)y-tip_local[1])/radius);
      double score = spacefield_(x, y) - 0.25*dist2; // slight preference for result near the centroid
      if (score > best_score)
      {
        best_score = score;
        best_x = x;
        best_y = y;
      }
    }
  }
  if (best_score > 0.0) 
  {
    tip[0] = ((double)best_x+0.5)*voxel_width_;
    tip[1] = ((double)best_y+0.5)*voxel_width_;
    return true;
  }
  return false;
}

// extract the ray cloud canopy to a height field, then call the heightfield based forest extraction
ray::ForestStructure Forest::extract(const std::string &cloud_name_stub, Mesh &mesh, const std::vector<std::pair<Eigen::Vector3d, double> > &trunks, double voxel_width) 
{
  trunks_ = trunks;
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name_stub + ".ply", info))
    return ray::ForestStructure();
  min_bounds_ = info.ends_bound.min_bound_;
  max_bounds_ = info.ends_bound.max_bound_;

  std::cout << "voxel width: " << voxel_width << " m" << std::endl;

  double width = (max_bounds_[0] - min_bounds_[0])/voxel_width;
  double length = (max_bounds_[1] - min_bounds_[1])/voxel_width;
  Eigen::Vector2i grid_dims(ceil(width), ceil(length));
  std::cout << "dims for heightfield: " << grid_dims.transpose() << std::endl;
  Eigen::ArrayXXd highs = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], -1e10);

  auto fillHeightField = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends, std::vector<double> &, std::vector<ray::RGBA> &colours)
  {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha == 0)
        continue;
      Eigen::Vector3d pos = (ends[i] - min_bounds_)/voxel_width;
      double &h = highs((int)pos[0], (int)pos[1]);
      h = std::max(h, ends[i][2]);
    }
  };
  if (!ray::Cloud::read(cloud_name_stub + ".ply", fillHeightField))
    return ray::ForestStructure();

  Eigen::ArrayXXd lows;
  if (mesh.vertices().empty())
    lows = Eigen::ArrayXXd::Constant(highs.rows(), highs.cols(), min_bounds_[2]);
  else 
    mesh.toHeightField(lows, min_bounds_, max_bounds_, voxel_width);
  if (lows.rows() != highs.rows() || lows.cols() != highs.cols())
    std::cerr << "error: arrays are different widths " << lows.rows() << "!=" << highs.rows() << " or " << lows.cols() << "!=" << highs.cols() << std::endl;

  // generate grid
  Occupancy2D grid2D;
  if (!grid2D.load(cloud_name_stub + "_occupied.dat"))
  {
    grid2D.init(min_bounds_, max_bounds_, voxel_width);
    // walk the rays to fill densities
    grid2D.fillDensities(cloud_name_stub + ".ply", lows, 1.0, 1.5);
    grid2D.save(cloud_name_stub + "_occupied.dat");
  }
  if (grid2D.dims_[0] != lows.rows() || grid2D.dims_[1] != lows.cols())
    std::cerr << "error: arrays are different widths " << lows.rows() << "!=" << grid2D.dims_[0] << " or " << lows.cols() << "!=" << grid2D.dims_[1] << std::endl;
  Eigen::ArrayXXd space(grid2D.dims_[0], grid2D.dims_[1]);
  for (int i = 0; i<space.rows(); i++)
  {
    for (int j = 0; j<space.cols(); j++)
      space(i,j) = grid2D.pixel(Eigen::Vector3i(i, j, 0)).density();
  }

  return extract(highs, lows, space, voxel_width, cloud_name_stub);
}

void Forest::addTrunkHeights()
{
  for (int c = 0; c<(int)trunks_.size(); c++) // if there are known trunks, then include them...
  {
    auto &trunk = trunks_[c];
    Eigen::Vector3d posr = (trunk.first - min_bounds_)/voxel_width_;
    Eigen::Vector3i pos = posr.cast<int>();
    if (pos[0] < 0 || pos[0] >= heightfield_.rows() || pos[1] < 0 || pos[1] >= heightfield_.cols())
    {
      continue;
    }
    double radius = 10.0 * trunk.second / voxel_width_;
    double height = 80.0 * trunk.second;
    int rad = (int)std::ceil(radius);
    for (int x = std::max(0, pos[0]-rad); x <= std::min(pos[0] + rad, (int)heightfield_.rows()-1); x++)
    {
      for (int y = std::max(0, pos[1]-rad); y <= std::min(pos[1] + rad, (int)heightfield_.cols()-1); y++)
      {
        Eigen::Vector2d dif((double)x + 0.5 - posr[0], (double)y + 0.5 - posr[1]);
        dif /= radius;
        double r = dif.squaredNorm();
        double h = height * (1.0 - r);
        if (h > 0.0)
        {
          if (heightfield_(x, y) != -1e10)
            heightfield_(x, y) += h;
        }
      }    
    }
  }
}

void Forest::smoothHeightfield()
{
  Eigen::ArrayXXd smooth_heights = heightfield_;
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
      double &h = heightfield_(x, y);
      if (h == -1e10)
        continue;
      double mean = h;
      double count = 1;
      for (int xx = std::max(0, x-1); xx <= std::min(x+1, (int)heightfield_.rows()-1); xx++)
      {
        for (int yy = std::max(0, y-1); yy <= std::min(y+1, (int)heightfield_.cols()-1); yy++)
        {
          double &h2 = heightfield_(xx, yy);
          if (h2 != -1e10)
          {
            mean += h2;
            count++;
          }
        }        
      }
      smooth_heights(x,y) = mean / count;
    }
  }
  heightfield_ = smooth_heights;
}

ray::ForestStructure Forest::extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, const Eigen::ArrayXXd &space, double voxel_width, const std::string &cloud_name_stub)
{
  voxel_width_ = voxel_width;
  heightfield_ = highs;
  lowfield_ = lows;
  spacefield_ = space;
  drawHeightField(cloud_name_stub + "_highfield.png", heightfield_);
  drawHeightField(cloud_name_stub + "_lowfield.png", lowfield_);
   
  
  ray::ForestStructure forest;
  int num_spaces = 0;
  const std::vector<std::string> attributes = {"tree_radius", "height", "trunk_identified"};
  int tree_radius_id = 0;
  int height_id = 1;
  int trunk_identified_id = 2;

  if (agglomerate_)
  {
    int count = 0;
    std::vector<Eigen::Vector3d> points;
    for (int x = 0; x < heightfield_.rows(); x++)
    {
      for (int y = 0; y < heightfield_.cols(); y++)
      {
        double &h = heightfield_(x, y);
        double l = lowfield_(x, y);
        if (h < l+undercroft_height)
        {
          h = -1e10;
          count++;
        }
        if (h > -1e10 && l > -1e10 && h>=l)
          points.push_back(Eigen::Vector3d(voxel_width_*((double)x + 0.5), voxel_width_*((double)y + 0.5), h-l)); // get heightfield relative to ground
      }
    }
    std::cout << "undercroft at height " << undercroft_height << " removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;

    Terrain terrain;
    terrain.growDownwards(points, wrap_gradient);
    Mesh &mesh = terrain.mesh();

    std::cout << "num points " << mesh.vertices().size() << std::endl;
    mesh.reduce();
    std::cout << "num verts: " << mesh.vertices().size() << std::endl;


    // 2. cluster according to radius based on height of points
    std::vector<Cluster> point_clusters;
    agglomerate(mesh.vertices(), mesh.indexList(), min_diameter_per_height_, max_diameter_per_height_, point_clusters);
    std::cout << "number found: " << point_clusters.size() << std::endl;
    std::vector<Eigen::Vector3d> &verts = mesh.vertices();

    if (verbose)
    {
      renderAgglomeration(point_clusters, verts, cloud_name_stub);
      drawAgglomeration(point_clusters, verts, cloud_name_stub);
    }
    
    for (auto &cluster: point_clusters)
    {
      Eigen::Vector3d tip;
      int trunk_id = cluster.trunk_id;    
      if (findSpace(cluster, verts, tip))
      {
        ray::TreeStructure tree;
        ray::TreeStructure::Segment result;
        tree.attributes() = attributes;

        result.tip = min_bounds_ + tip;
        result.tip[2] = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
        result.attributes[height_id] = tip[2]; 
        result.attributes[trunk_identified_id] = 1; 
        if (trunk_id >= 0)
          result.radius = trunks_[trunk_id].second;
        else 
        {
          result.radius = result.attributes[height_id] / approx_height_per_radius_;
          result.attributes[trunk_identified_id] = 0;
        }
        tree.segments().push_back(result);
        forest.trees.push_back(tree);
      }
      else
      {
        num_spaces++;
      }
    }  
  }
  else
  {
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
    original_heightfield_ = heightfield_;
    addTrunkHeights();
    drawHeightField(cloud_name_stub + "_trunkhighfield.png", heightfield_);
    for (int i = 0; i<smooth_iterations_; i++)
      smoothHeightfield();
    drawHeightField(cloud_name_stub + "_smoothhighfield.png", heightfield_);

    std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;
    std::vector<TreeNode> trees;
    std::set<int> heads;
    hierarchicalWatershed(trees, heads);

    std::cout << "number of raw candidates: " << trees.size() << " number largest size: " << heads.size() << std::endl;

    // calculate the area of pixels occupied by each index
    for (int x = 0; x < indexfield_.rows(); x++)
    {
      for (int y = 0; y < indexfield_.cols(); y++)
      {
        int ind = indexfield_(x, y);
        if (ind == -1)
          continue;
        while (trees[ind].attaches_to != -1)
          ind = trees[ind].attaches_to;
        trees[ind].area++;
      }
    }

    drawFinalSegmentation(cloud_name_stub, trees);
    renderWatershed(cloud_name_stub, trees, heads);

    for (auto &ind: heads)
    {
      if (trees[ind].area < min_area_)
        continue;
      Eigen::Vector3d tip;
      int trunk_id = trees[ind].trunk_id;
      if (findSpace2(trees[ind], tip))
      {
        ray::TreeStructure tree;
        tree.attributes() = attributes;
        ray::TreeStructure::Segment result;
        result.tip = min_bounds_ + tip;
        result.tip[2] = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
        result.attributes[height_id] = tip[2]; 
        int num_pixels = trees[ind].area;
        result.attributes[tree_radius_id] = std::sqrt(((double)num_pixels * voxel_width_*voxel_width_)/kPi); // get from num pixels
        result.attributes[trunk_identified_id] = 1; 
        if (trunk_id >= 0)
          result.radius = trunks_[trunk_id].second;
        else 
        {
          result.radius = result.attributes[height_id] / approx_height_per_radius_;
          result.attributes[trunk_identified_id] = 0; 
        }
        tree.segments().push_back(result);
        forest.trees.push_back(tree);
      }
      else
      {
        num_spaces++;
      }
    }    
  }
  std::cout << "number of disallowed trees: " << num_spaces << " / " << forest.trees.size() << std::endl;

  std::sort(forest.trees.begin(), forest.trees.end(), [](const ray::TreeStructure &a, const ray::TreeStructure &b){ return a.segments()[0].attributes[0] > b.segments()[0].attributes[0]; });

  return forest;
}


}
