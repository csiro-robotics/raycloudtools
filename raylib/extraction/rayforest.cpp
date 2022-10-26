// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"
#include "../rayconvexhull.h"
#include "../raycuboid.h"
#include "../rayforestgen.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "rayterrain.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace ray
{
/// find space in the forest in which a trunk could reside
bool Forest::findSpace(const TreeNode &node, Eigen::Vector3d &tip)
{
  bool calculate = true;
  // set the tip as the mid points of the node's bounds
  Eigen::Vector2d mid = voxel_width_ * (node.min_bound + node.max_bound).cast<double>() / 2.0;
  tip[0] = mid[0];
  tip[1] = mid[1];

  // if this node is associated with a trunk, then use the trunk location, not the centroid
  if (node.trunk_id >= 0)
  {
    tip = trunks_[node.trunk_id].first - min_bounds_;
    calculate = false;
  }

  const Eigen::Vector3d tip_local = tip / voxel_width_;
  tip[2] = node.peak[2] - lowfield_(static_cast<int>(tip_local[0]), static_cast<int>(tip_local[1]));
  if (!calculate)
  {
    return true;
  }
  // use a fixed cone angle for the downwards search
  const double search_down_gradient = 0.2;
  const double radius = tip_local[2] * search_down_gradient;

  // now find the closest bit of space to put the tree in:
  const int min_x = std::max(0, static_cast<int>(tip_local[0] - radius));
  const int max_x = std::min(static_cast<int>(spacefield_.rows()) - 1, static_cast<int>(tip_local[0] + radius));
  const int min_y = std::max(0, static_cast<int>(tip_local[1] - radius));
  const int max_y = std::min(static_cast<int>(spacefield_.cols()) - 1, static_cast<int>(tip_local[1] + radius));
  double best_score = std::numeric_limits<double>::lowest();
  int best_x = -1;
  int best_y = -1;
  // for each cell with in the calculated bounds
  for (int x = min_x; x <= max_x; x++)
  {
    for (int y = min_y; y <= max_y; y++)
    {
      const double dist2 =
        sqr((static_cast<double>(x) - tip_local[0]) / radius) + sqr((static_cast<double>(y) - tip_local[1]) / radius);
      const double score = spacefield_(x, y) - 0.25 * dist2;  // slight preference for result near the centroid
      if (score > best_score)
      {
        best_score = score;
        best_x = x;
        best_y = y;
      }
    }
  }
  // choose the point that has space (score > 0) and that has the best score
  if (best_score > 0.0)
  {
    tip[0] = (static_cast<double>(best_x) + 0.5) * voxel_width_;
    tip[1] = (static_cast<double>(best_y) + 0.5) * voxel_width_;
    return true;
  }
  return false;
}

// extract the ray cloud canopy to a height field, then call the heightfield based forest extraction
ray::ForestStructure Forest::extract(const std::string &cloud_name_stub, Mesh &mesh,
                                     const std::vector<std::pair<Eigen::Vector3d, double>> &trunks, double voxel_width)
{
  trunks_ = trunks;
  // firstly, get the bounds of the ray cloud
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name_stub + ".ply", info))
  {
    return ray::ForestStructure();
  }
  min_bounds_ = info.ends_bound.min_bound_;
  max_bounds_ = info.ends_bound.max_bound_;

  // then we need to generate some height fields, these are 2D arrays
  const double width = (max_bounds_[0] - min_bounds_[0]) / voxel_width;
  const double length = (max_bounds_[1] - min_bounds_[1]) / voxel_width;
  const Eigen::Vector2i grid_dims(ceil(width), ceil(length));
  std::cout << "dims for heightfield: " << grid_dims.transpose() << std::endl;
  Eigen::ArrayXXd highs = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], std::numeric_limits<double>::lowest());

  // fill in the highest points on the input cloud
  auto fillHeightField = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                             std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha == 0)
      {
        continue;
      }
      const Eigen::Vector3d pos = (ends[i] - min_bounds_) / voxel_width;
      double &h = highs(static_cast<int>(pos[0]), static_cast<int>(pos[1]));
      h = std::max(h, ends[i][2]);
    }
  };
  if (!ray::Cloud::read(cloud_name_stub + ".ply", fillHeightField))
  {
    return ray::ForestStructure();
  }

  // next fill in the lowest points using the supplied ground mesh
  Eigen::ArrayXXd lows;
  if (mesh.vertices().empty())
  {
    lows = Eigen::ArrayXXd::Constant(highs.rows(), highs.cols(), min_bounds_[2]);
  }
  else
  {
    mesh.toHeightField(lows, min_bounds_, max_bounds_, voxel_width);
  }
  if (lows.rows() != highs.rows() || lows.cols() != highs.cols())
  {
    std::cerr << "error: arrays are different widths " << lows.rows() << "!=" << highs.rows() << " or " << lows.cols()
              << "!=" << highs.cols() << std::endl;
  }

  // generate a 2D grid in order to fill in the 'space field' a 2D array of free space (where the rays are)
  OccupancyGrid2D grid2D;
  if (!grid2D.load(cloud_name_stub + "_occupied.dat"))
  {
    grid2D.init(min_bounds_, max_bounds_, voxel_width);
    // walk the rays to fill densities based on walking the rays through the grid
    grid2D.fillDensities(cloud_name_stub + ".ply", lows, 1.0, 1.5);
    grid2D.save(cloud_name_stub + "_occupied.dat");
  }
  if (grid2D.dims()[0] != lows.rows() || grid2D.dims()[1] != lows.cols())
  {
    std::cerr << "error: arrays are different widths " << lows.rows() << "!=" << grid2D.dims()[0] << " or "
              << lows.cols() << "!=" << grid2D.dims()[1] << std::endl;
  }
  // move the grid into a 2D 'space' array
  Eigen::ArrayXXd space(grid2D.dims()[0], grid2D.dims()[1]);
  for (int i = 0; i < space.rows(); i++)
  {
    for (int j = 0; j < space.cols(); j++)
    {
      space(i, j) = grid2D.pixel(Eigen::Vector3i(i, j, 0)).density();
    }
  }

  return extract(highs, lows, space, voxel_width, cloud_name_stub);
}

/// include any previously observed trunks into the height field as bumps (paraboloids)
/// this is a soft hint for where the trees should be found
void Forest::addTrunkHeights()
{
  for (int c = 0; c < static_cast<int>(trunks_.size()); c++)  // if there are known trunks, then include them...
  {
    auto &trunk = trunks_[c];
    const Eigen::Vector3d posr = (trunk.first - min_bounds_) / voxel_width_;
    const Eigen::Vector3i pos = posr.cast<int>();
    if (pos[0] < 0 || pos[0] >= heightfield_.rows() || pos[1] < 0 || pos[1] >= heightfield_.cols())
    {
      continue;
    }
    // an approximate radius around the trunk
    const double radius = 10.0 * trunk.second / voxel_width_;
    // define a steepness for the bump due to the trunk
    const double height = 80.0 * trunk.second;
    const int rad = static_cast<int>(std::ceil(radius));
    for (int x = std::max(0, pos[0] - rad); x <= std::min(pos[0] + rad, static_cast<int>(heightfield_.rows()) - 1); x++)
    {
      for (int y = std::max(0, pos[1] - rad); y <= std::min(pos[1] + rad, static_cast<int>(heightfield_.cols()) - 1);
           y++)
      {
        Eigen::Vector2d dif(static_cast<double>(x) + 0.5 - posr[0], static_cast<double>(y) + 0.5 - posr[1]);
        dif /= radius;
        const double r = dif.squaredNorm();
        const double h = height * (1.0 - r);  // this is the paraboloid height
        if (h > 0.0)
        {
          if (heightfield_(x, y) != std::numeric_limits<double>::lowest())
          {
            heightfield_(x, y) += h;
          }
        }
      }
    }
  }
}

/// Smooth the height field to remove noise that can interfere with the signal of tree crowns
void Forest::smoothHeightfield()
{
  Eigen::ArrayXXd smooth_heights = heightfield_;
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
      double &h = heightfield_(x, y);
      if (h == std::numeric_limits<double>::lowest())
      {
        continue;
      }
      double mean = h;
      double count = 1;
      // use a mean of the valid heights in the Moore neighbourhood of each cell
      for (int xx = std::max(0, x - 1); xx <= std::min(x + 1, static_cast<int>(heightfield_.rows()) - 1); xx++)
      {
        for (int yy = std::max(0, y - 1); yy <= std::min(y + 1, static_cast<int>(heightfield_.cols()) - 1); yy++)
        {
          const double &h2 = heightfield_(xx, yy);
          if (h2 != std::numeric_limits<double>::lowest())
          {
            mean += h2;
            count++;
          }
        }
      }
      smooth_heights(x, y) = mean / count;
    }
  }
  heightfield_ = smooth_heights;
}

/// extract tree locations from a set of three 2D arrays, a height field (the canopy) a low field (the ground)
/// and a space field (the free space)
ray::ForestStructure Forest::extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows,
                                     const Eigen::ArrayXXd &space, double voxel_width,
                                     const std::string &cloud_name_stub)
{
  voxel_width_ = voxel_width;
  heightfield_ = highs;
  lowfield_ = lows;
  spacefield_ = space;
  // useful debug drawing
  drawHeightField(cloud_name_stub + "_highfield.png", heightfield_);
  drawHeightField(cloud_name_stub + "_lowfield.png", lowfield_);

  // the output is a fores structure
  ray::ForestStructure forest;
  int num_spaces = 0;
  // and we include four user-defined attributes
  const std::vector<std::string> tree_attributes = { "tree_radius", "height", "trunk_identified" };
  const int tree_radius_id = 0;
  const int height_id = 1;
  const int trunk_identified_id = 2;

  // the indexfield assigns a unique index to each cluster (each tree) as we grow using the watershed algorithm
  indexfield_ = Eigen::ArrayXXi::Constant(heightfield_.rows(), heightfield_.cols(), -1);

  // ignore the undercroft
  int count = 0;
  for (int x = 0; x < heightfield_.rows(); x++)
  {
    for (int y = 0; y < heightfield_.cols(); y++)
    {
      if (heightfield_(x, y) < lowfield_(x, y) + undercroft_height)
      {
        heightfield_(x, y) = std::numeric_limits<double>::lowest();
        count++;
      }
    }
  }

  original_heightfield_ = heightfield_;
  addTrunkHeights();
  drawHeightField(cloud_name_stub + "_trunkhighfield.png", heightfield_);
  for (int i = 0; i < smooth_iterations_; i++)
  {
    smoothHeightfield();
  }
  drawHeightField(cloud_name_stub + "_smoothhighfield.png", heightfield_);

  std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows() * heightfield_.cols() << std::endl;
  std::vector<TreeNode> trees;
  std::set<int> heads;

  // this is the main segmentation algorithm, it generates the trees vector
  hierarchicalWatershed(trees, heads);

  std::cout << "number of raw candidates: " << trees.size() << " number largest size: " << heads.size() << std::endl;

  // calculate the area of pixels occupied by each index
  for (int x = 0; x < indexfield_.rows(); x++)
  {
    for (int y = 0; y < indexfield_.cols(); y++)
    {
      int ind = indexfield_(x, y);
      if (ind == -1)
      {
        continue;
      }
      while (trees[ind].attaches_to != -1)
      {
        ind = trees[ind].attaches_to;
      }
      trees[ind].area++;
    }
  }

  drawFinalSegmentation(cloud_name_stub, trees);
  renderWatershed(cloud_name_stub, trees, heads);

  // now we generate the actual output 'forest' structure
  for (auto &ind : heads)  // for each tree
  {
    if (trees[ind].area < min_area_)  // too small to count as a tree
    {
      continue;
    }
    Eigen::Vector3d tip;
    const int trunk_id = trees[ind].trunk_id;
    if (findSpace(trees[ind], tip))  // if this tree actually has space to exist
    {
      ray::TreeStructure tree;
      tree.treeAttributeNames() = tree_attributes;
      tree.attributeNames().push_back("section_id");
      tree.treeAttributes().resize(tree_attributes.size());
      ray::TreeStructure::Segment result;
      // locate the tree
      result.tip = min_bounds_ + tip;
      result.tip[2] = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
      // set its height
      tree.treeAttributes()[height_id] = tip[2];
      // estimate the tree (crown) radius
      const int num_pixels = trees[ind].area;
      tree.treeAttributes()[tree_radius_id] =
        std::sqrt((static_cast<double>(num_pixels) * voxel_width_ * voxel_width_) / kPi);  // get from num pixels
      tree.treeAttributes()[trunk_identified_id] = 1;
      // assign its unique section id
      result.attributes.push_back(ind);
      // if the tree had an identified trunk then use this radius estimate
      if (trunk_id >= 0)
      {
        result.radius = trunks_[trunk_id].second;
      }
      else  // otherwise, estimate trunk radius crudely from the tree height
      {
        result.radius = tree.treeAttributes()[height_id] / approx_height_per_radius_;
        tree.treeAttributes()[trunk_identified_id] = 0;
      }
      tree.segments().push_back(result);
      forest.trees.push_back(tree);
    }
    else
    {
      num_spaces++;
    }
  }
  std::cout << "number of disallowed trees: " << num_spaces << " / " << forest.trees.size() << std::endl;

  // sort trees by their radius
  std::sort(forest.trees.begin(), forest.trees.end(),
            [&tree_radius_id](const ray::TreeStructure &a, const ray::TreeStructure &b) {
              return a.treeAttributes()[tree_radius_id] > b.treeAttributes()[tree_radius_id];
            });

  return forest;
}


}  // namespace ray
