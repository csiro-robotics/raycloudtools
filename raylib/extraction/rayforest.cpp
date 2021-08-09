// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"
#include "../rayconvexhull.h"
#include "../raymesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <nabo/nabo.h>

namespace ray
{
void agglomerate(const std::vector<Eigen::Vector3d> &points, double min_diameter, double max_diameter, std::vector< std::vector<int> > &point_clusters)
{
  // 1. get nearest neighbours for each point
  const int search_size = std::min(8, (int)points.size()-1);
  Eigen::MatrixXd points_p(2, points.size());
  for (unsigned int i = 0; i < points.size(); i++) 
    points_p.col(i) = Eigen::Vector2d(points[i][0], points[i][1]);
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 2);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, min_diameter);

  struct Nd
  {
    Nd(int id1, int id2, double dist2) : id1(id1), id2(id2), dist2(dist2) {}
    int id1, id2;
    double dist2;
  };
  std::vector<Nd> nds;
  for (size_t i = 0; i<points.size(); i++)
  {
    for (int j = 0; j<search_size && indices(j, i) > -1; j++)
    {
      nds.push_back(Nd((int)i, indices(j, i), dists2(j, i)));
    }
  }
  std::sort(nds.begin(), nds.end(), [](const Nd &nd1, const Nd &nd2){ return nd1.dist2 < nd2.dist2; });
  struct Cluster
  {
    Eigen::Vector3d min_bound, max_bound;
    std::vector<int> ids;
    bool active;
  };
  std::vector<Cluster> clusters(points.size());
  std::vector<int> cluster_ids(points.size());
  std::vector<bool> visited(points.size());
  for (size_t i = 0; i<points.size(); i++)
  {
    clusters[i].min_bound = clusters[i].max_bound = points[i];
    clusters[i].active = true;
    clusters[i].ids.push_back((int)i);
    cluster_ids[i] = (int)i;
    visited[i] = false;
  }

  // 2. for each node in turn, from smallest to highest distance, agglomerate
  for (auto &node: nds)
  {
    if (cluster_ids[node.id1] == cluster_ids[node.id2]) // already part of same cluster
      continue; 
    int cl1 = cluster_ids[node.id1];
    int cl2 = cluster_ids[node.id2];
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = std::max(dims[0], dims[1]);
    if (diam < max_diameter) // then merge
    {
      int first = std::min(cl1, cl2);
      int last = std::max(cl1, cl2);
      clusters[first].min_bound = minb;
      clusters[first].max_bound = maxb;
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
      point_clusters.push_back(cluster.ids);
    }
  }
}



// extract the ray cloud canopy to a height field, then call the heightfield based forest extraction
void Forest::extract(const Cloud &cloud, Mesh &mesh, double voxel_width) // = 4.0 * cloud.estimatePointSpacing();
{
  cloud.calcBounds(&min_bounds_, &max_bounds_);

  double width = (max_bounds_[0] - min_bounds_[0])/voxel_width;
  double length = (max_bounds_[1] - min_bounds_[1])/voxel_width;
  Eigen::Vector2i grid_dims(ceil(width), ceil(length));
  std::cout << "dims for heightfield: " << grid_dims.transpose() << std::endl;
  Eigen::ArrayXXd highs = Eigen::ArrayXXd::Constant(grid_dims[0], grid_dims[1], -1e10);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    const Eigen::Vector3d &p = cloud.ends[i];
    Eigen::Vector3d pos = (p - min_bounds_)/voxel_width;
    double &h = highs((int)pos[0], (int)pos[1]);
    h = std::max(h, p[2]);
  }
  Eigen::ArrayXXd lows;
  mesh.toHeightField(lows, min_bounds_, max_bounds_, voxel_width);
  extract(highs, lows, voxel_width);
}

void Forest::extract(const Eigen::ArrayXXd &highs, const Eigen::ArrayXXd &lows, double voxel_width)
{
  voxel_width_ = voxel_width;
  heightfield_ = highs;
  lowfield_ = lows;
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
      if (h > -1e10 && l > -1e10)
      {
        points.push_back(Eigen::Vector3d(voxel_width_*(double)x, voxel_width_*(double)y, h-l)); // get heightfield relative to ground
      }
    }
  }
  std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;
  drawHeightField("highfield.png", heightfield_);
  drawHeightField("lowfield.png", lowfield_);

  const double curvature = 0.1;
  const double crown_radius_per_height = 0.3;
  // 1. get top-down mesh of canopy
  ConvexHull hull(points);
  hull.growDownwards(curvature);
  Mesh &mesh = hull.mesh();
  std::vector<Eigen::Vector3d> &verts = mesh.vertices();
  //std::vector<Eigen::Vector3i> &inds = mesh.indexList();

  // 2. cluster according to radius based on height of points
  std::vector< std::vector<int> > point_clusters;
  double min_diameter = 0.5; 
  double max_diameter = 5.0;
  agglomerate(verts, min_diameter, max_diameter, point_clusters);

  // 3. for each cluster, get a mean centre...
  for (auto &cluster: point_clusters)
  {
    Eigen::Vector3d weighted_sum(0,0,0);
    double weight = 0.0;
    for (auto &i: cluster)
    {
      weighted_sum += verts[i][2] * verts[i];
      weight += verts[i][2];
    }
    Result tree;
    Eigen::Vector3d tip = weighted_sum / weight;
    tree.ground_height = lowfield_(int(tip[0] / voxel_width_), int(tip[1] / voxel_width_));
    tree.tree_tip = tip + Eigen::Vector3d(0,0,tree.ground_height);
    tree.radius = tip[2] * crown_radius_per_height;
    tree.curvature = 1.0; // TODO fix
    results_.push_back(tree);
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
    Eigen::Vector3d base = result.tree_tip * voxel_width_;
    base[2] = result.ground_height;
    const double tree_radius_to_trunk_radius = 1.0/20.0; // TODO: temporary until we have a better parameter choice
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << result.radius*tree_radius_to_trunk_radius << std::endl;
  }
  return true;
}
}
