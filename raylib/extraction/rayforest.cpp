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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#define PARABOLOID

namespace ray
{

bool Forest::findSpace(const Cluster &cluster, const std::vector<Eigen::Vector3d> &points, Eigen::Vector3d &tip)
{
  if (cluster.ids.empty())
  {
    tip = cluster.max_bound;
    return true; // TODO: do we believe truks absolutely? or should we double check against the free space?
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
  const double gradient = 0.1;
  double radius = tip_local[2]*gradient;

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
  if (best_x != -1) 
  {
    tip[0] = ((double)best_x+0.5)*voxel_width_;
    tip[1] = ((double)best_y+0.5)*voxel_width_;
    return true;
  }
  return false;
}

void Forest::agglomerate(const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> index_list, double min_diameter_per_height, double max_diameter_per_height, std::vector<Cluster> &point_clusters)
{
  struct Nd
  {
    Nd(int id1, int id2, double dist2) : id1(id1), id2(id2), dist2(dist2) {}
    int id1, id2;
    double dist2;
  };
  std::vector<Nd> nds;
  for (auto &ind: index_list) // doubles up on edges, which slows down the sort, but is safely discarded in later code
  {
    for (int i = 0; i<3; i++)
    {
      int id1 = ind[i];
      int id2 = ind[(i+1)%3];
      Eigen::Vector3d diff = points[id1]-points[id2];
      #if defined PARABOLOID // paraboloid allows points sweeping down the side of trees, so shouldn't penalise vertical distance
      diff[2] = 0.0;
      #endif
      nds.push_back(Nd(std::min(id1, id2), std::max(id1, id2), diff.squaredNorm()));
    }
  }
  std::sort(nds.begin(), nds.end(), [](const Nd &nd1, const Nd &nd2){ return nd1.dist2 < nd2.dist2; });

  std::vector<Cluster> clusters(points.size());
  std::vector<int> cluster_ids(points.size());
  std::vector<bool> visited(points.size());
  for (size_t i = 0; i<points.size(); i++)
  {
    clusters[i].min_bound = clusters[i].max_bound = points[i];
    clusters[i].active = true;
    clusters[i].ids.push_back((int)i);
    clusters[i].trunk_id = -1;
    cluster_ids[i] = (int)i;
    visited[i] = false;
  }
  int outside_count = 0;
  for (int c = 0; c<(int)trunks_.size(); c++) // if there are known trunks, then include them...
  {
    auto &trunk = trunks_[c];
    double min_dist2 = 1e10f;
    int closest_i = -1;
    for (int i = 0; i<(int)points.size(); i++)
    {
      Eigen::Vector3d dif = points[i] - (trunk.first - min_bounds_);
      dif[2] = 0.0;
      double dist2 = dif.squaredNorm();
      if (dist2 < min_dist2)
      {
        min_dist2 = dist2;
        closest_i = i;
      }
    }
    const double gradient = 0.1;
    if (closest_i >= 0 && std::sqrt(min_dist2) < gradient*points[closest_i][2])
      clusters[cluster_ids[closest_i]].trunk_id = c;
    else // this trunk has no peak points above it, so is a lower tree. We still need to calculate its height
    {
      // this uses heightfield
      int rad = 3;
      Eigen::Vector3i index = ((trunk.first-min_bounds_)/voxel_width_).cast<int>(); 
      int minx = std::max(0, index[0]-rad);
      int maxx = std::min(index[0]+rad, (int)heightfield_.rows()-1);
      int miny = std::max(0, index[1]-rad);
      int maxy = std::min(index[1]+rad, (int)heightfield_.cols()-1);
      double maxheight = -1e10;
      for (int x = minx; x<=maxx; x++)
        for (int y = miny; y<=maxy; y++)
          maxheight = std::max(maxheight, heightfield_(x, y));
      std::cout << "trunk " << ++outside_count << "/" << trunks_.size() << " has no peaks above it, canopy found at " << maxheight << " above trunk" << std::endl;
      // now how to we relay this information??
      Eigen::Vector3d tip = trunk.first - min_bounds_;
      tip[2] = maxheight;
      Cluster new_cluster;
      new_cluster.min_bound = new_cluster.max_bound = tip;
      new_cluster.trunk_id = c;
      clusters.push_back(new_cluster);
    }
  }

  // 2. for each node in turn, from smallest to highest distance, agglomerate
  for (auto &node: nds)
  {
    int cl1 = cluster_ids[node.id1];
    int cl2 = cluster_ids[node.id2];
    if (cl1 == cl2) // already part of same cluster
      continue; 
    if (clusters[cl1].trunk_id >=0 && clusters[cl2].trunk_id >=0 && clusters[cl1].trunk_id != clusters[cl2].trunk_id) // don't merge from different trunks
    {
      std::cout << "clusters " << cl1 << " and " << cl2 << "cannot merge between two different trunks " << clusters[cl1].trunk_id << " and " << clusters[cl2].trunk_id << std::endl;
      continue;
    } 
    if (clusters[cl1].trunk_id >=0 && clusters[cl1].trunk_id == clusters[cl2].trunk_id)
      std::cout << "weird, two different clusters have the same trunk id " << cl1 << ", " << cl2 << " have " << clusters[cl1].trunk_id << " and " << clusters[cl2].trunk_id << std::endl;
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = std::max(dims[0], dims[1]);
    double mean_height = (minb[2] + maxb[2])/2.0; 

    double dist2 = node.dist2 / (points[node.id1][2]*points[node.id2][2]);
    if (dist2 > sqr(min_diameter_per_height))
    {
      // about to keep as isolated cluster. So check if there is space underneath, and only isolate if there is
      //#define MERGE_IF_NO_SPACE
      #if defined MERGE_IF_NO_SPACE 
      Eigen::Vector3d tip;
      bool space1 = findSpace(clusters[cl1], points, tip);
      bool space2 = findSpace(clusters[cl2], points, tip);
      if (space1 && space2) // space is found for each so they're separate, the end
      #endif
        continue; 
    }
    if (diam < max_diameter_per_height * mean_height) // then merge
    {
      int first = std::min(cl1, cl2);
      int last = std::max(cl1, cl2);
      clusters[first].min_bound = minb;
      clusters[first].max_bound = maxb;
      if (clusters[last].trunk_id >= 0)
      {
        if (clusters[first].trunk_id >= 0)
          std::cout << "error, this should never happen " << clusters[first].trunk_id << ", " << clusters[last].trunk_id << std::endl;
        clusters[first].trunk_id = clusters[last].trunk_id;
      }
      clusters[first].ids.insert(clusters[first].ids.begin(), clusters[last].ids.begin(), clusters[last].ids.end());
      for (auto &id: clusters[last].ids)
        cluster_ids[id] = first;
      clusters[last].active = false;
    }
  }
  for (auto &cluster: clusters)
  {
    if (cluster.active)
    {
      point_clusters.push_back(cluster);
    }
  }
}



// extract the ray cloud canopy to a height field, then call the heightfield based forest extraction
bool Forest::extract(const std::string &cloud_name, Mesh &mesh, const std::vector<std::pair<Eigen::Vector3d, double> > &trunks) 
{
  trunks_ = trunks;
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
    return false;
  min_bounds_ = info.ends_bound.min_bound_;
  max_bounds_ = info.ends_bound.max_bound_;
  double voxel_width = 0.25; // 6.0 * Cloud::estimatePointSpacing(cloud_name, info.ends_bound, info.num_bounded);
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

  if (!ray::Cloud::read(cloud_name, fillHeightField))
    return false;

  Eigen::ArrayXXd lows;
  mesh.toHeightField(lows, min_bounds_, max_bounds_, voxel_width);
  if (lows.rows() != highs.rows() || lows.cols() != highs.cols())
    std::cerr << "error: arrays are different widths" << std::endl;


    // generate grid
  Grid2D grid2D;
  if (!grid2D.load("occupied.dat"))
  {
    grid2D.init(info.ends_bound.min_bound_, info.ends_bound.max_bound_, voxel_width);
    // walk the rays to fill densities
    grid2D.fillDensities(cloud_name, lows, 1.0, 1.5);
    grid2D.save("occupied.dat");
  }
  if (grid2D.dims_[0] != lows.rows() || grid2D.dims_[1] != lows.cols())
    std::cerr << "error: arrays are different widths" << std::endl;
  Eigen::ArrayXXd space(grid2D.dims_[0], grid2D.dims_[1]);
  for (int i = 0; i<space.rows(); i++)
  {
    for (int j = 0; j<space.cols(); j++)
      space(i,j) = grid2D.pixel(Eigen::Vector3i(i, j, 0)).density();
  }

  extract(highs, lows, space, voxel_width);
  return true;
}

void Forest::extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, const Eigen::ArrayXXd &space, double voxel_width)
{
  voxel_width_ = voxel_width;
  heightfield_ = highs;
  lowfield_ = lows;
  spacefield_ = space;
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
  std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;
  drawHeightField("highfield.png", heightfield_);
  drawHeightField("lowfield.png", lowfield_);


#if defined PARABOLOID
  const double curvature_height = 4.0;
  const double max_diameter_per_height = 1.2; 
  const double min_diameter_per_height = 0.25;

  // now scale the points to maintain the curvature per height
  for (auto &point: points)
    point[2] = 0.5 * point[2]*point[2]; // squish
  // 1. get top-down mesh of canopy
  ConvexHull hull(points);
  hull.growDownwards(curvature_height);
  for (auto &point: points)
    point[2] = std::sqrt(2.0*point[2]); // unsquish
  Mesh &mesh = hull.mesh();
  mesh.vertices() = points;
#else
  const double max_diameter_per_height = 0.9;
  const double min_diameter_per_height = 0.1; 
  const double gradient = 1.0;

  Terrain terrain;
  terrain.growDownwards(points, gradient);
  Mesh &mesh = terrain.mesh();
#endif
  std::cout << "num points " << mesh.vertices().size() << std::endl;
  mesh.reduce();
  std::cout << "num verts: " << mesh.vertices().size() << std::endl;


  // 2. cluster according to radius based on height of points
  std::vector<Cluster> point_clusters;
  agglomerate(mesh.vertices(), mesh.indexList(), min_diameter_per_height, max_diameter_per_height, point_clusters);
  std::cout << "number found: " << point_clusters.size() << std::endl;
  std::vector<Eigen::Vector3d> &verts = mesh.vertices();
  const double height_per_radius = 50.0; // TODO: temporary until we have a better parameter choice

  if (verbose)
  {
    std::vector<Eigen::Vector3d> cloud_points;
    std::vector<double> times;
    std::vector<RGBA> colours;
    RGBA colour;
    colour.alpha = 255;
    for (auto &cluster: point_clusters)
    {
      colour.red = uint8_t(rand()%255);
      colour.green = uint8_t(rand()%255);
      colour.blue = uint8_t(rand()%255);
      for (auto &i: cluster.ids)
      {
        double ground_height = lowfield_(int(verts[i][0]/voxel_width_), int(verts[i][1]/voxel_width));
        Eigen::Vector3d pos = verts[i];
        pos += Eigen::Vector3d(min_bounds_[0], min_bounds_[1], ground_height + 0.05);
        cloud_points.push_back(pos);
        times.push_back(0.0);
        colours.push_back(colour);
      }
      Eigen::Vector3d tip;
      bool found = findSpace(cluster, verts, tip);
      double rad = tip[2] / height_per_radius;
      tip[2] = lowfield_(int(tip[0]/voxel_width_), int(tip[1]/voxel_width));
      tip[0] += min_bounds_[0];
      tip[1] += min_bounds_[1];
      if (!found)
        colour.red = colour.green = colour.blue = 0;
      double z_max = 2.0;
      if (cluster.trunk_id >= 0)
        z_max = 4.0;
      if (cluster.ids.empty())
      {
        colour.red = 255;
        colour.green = 0;
        colour.blue = 255;
      }
      for (double z = 0.0; z<z_max; z+=0.3)
      {
        for (double ang = 0.0; ang<2.0*kPi; ang += 0.3)
        {
          cloud_points.push_back(tip + Eigen::Vector3d(rad*std::sin(ang), rad*std::cos(ang), z));
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }

    // now add the space field:
    for (int i = 0; i<spacefield_.rows(); i++)
    {
      for (int j = 0; j<spacefield_.cols(); j++)
      {
        if (spacefield_(i,j) < 1.0)
        {
          double height = lowfield_(i, j) + 0.2;
          double x = min_bounds_[0] + (double)i*voxel_width_;
          double y = min_bounds_[1] + (double)j*voxel_width_;
          cloud_points.push_back(Eigen::Vector3d(x, y, height));
          times.push_back(0.0);
          colour.red = colour.green = colour.blue = (uint8_t)(255.0*spacefield_(i,j));
          colours.push_back(colour);
        }
      }
    }


    writePlyPointCloud("clusters.ply", cloud_points, times, colours);
  }

  // 3. for each cluster, get a mean centre...
  for (auto &cluster: point_clusters)
  {
    Eigen::Vector3d tip;
    if (findSpace(cluster, verts, tip))
    {
      Result tree;
      tree.height = tip[2];
      tree.base = min_bounds_ + tip;
      tree.base[2] = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
      if (cluster.trunk_id)
        tree.radius = trunks_[cluster.trunk_id].second;
      else 
        tree.radius = tip[2] / height_per_radius;
      results_.push_back(tree);
    }
  }

  drawTrees("result_trees.png", results_, (int)heightfield_.rows(), (int)heightfield_.cols());
}

bool Forest::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Forest extraction, tree base location list: x, y, z, radius" << std::endl;
  for (auto &result: results_)
  {
    ofs << result.base[0] << ", " << result.base[1] << ", " << result.base[2] << ", " << result.radius << ", " << result.height << std::endl;
  }
  return true;
}
}
