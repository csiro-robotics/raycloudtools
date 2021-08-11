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

namespace ray
{
void agglomerate(const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> index_list, double min_diameter_per_height, double max_diameter_per_height, std::vector< std::vector<int> > &point_clusters)
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
      diff[2] = 0.0;
      nds.push_back(Nd(std::min(id1, id2), std::max(id1, id2), diff.squaredNorm()));
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
    double dist2 = node.dist2 / (points[node.id1][2]*points[node.id2][2]);
    if (dist2 > sqr(min_diameter_per_height))
    {
 //     std::cout << "gap: " << std::sqrt(node.dist2) << " or " << (points[node.id1]-points[node.id2]).norm() << " at height: " << std::sqrt(points[node.id1][2]*points[node.id2][2]) << " > " << min_diameter_per_height << std::endl;
      continue;
    }
    int cl1 = cluster_ids[node.id1];
    int cl2 = cluster_ids[node.id2];
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = std::max(dims[0], dims[1]);
    double mean_height = (minb[2] + maxb[2])/2.0; 
    if (diam > 2.0*std::min(dims[0], dims[1]) && (clusters[cl1].ids.size() + clusters[cl2].ids.size()) > 4) // ignore merges that are too elongated. TODO: use eigenvalues eventually 
      continue;
    if (diam < max_diameter_per_height * mean_height) // then merge
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
bool Forest::extract(const std::string &cloud_name, Mesh &mesh) 
{
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
    return false;
  min_bounds_ = info.ends_bound.min_bound_;
  max_bounds_ = info.ends_bound.max_bound_;
  double voxel_width = 6.0 * Cloud::estimatePointSpacing(cloud_name, info.ends_bound, info.num_bounded);
  std::cout << "estimated voxel width: " << voxel_width << " m" << std::endl;

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
  extract(highs, lows, voxel_width);
  return true;
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
      if (h > -1e10 && l > -1e10 && h>=l)
        points.push_back(Eigen::Vector3d(voxel_width_*((double)x + 0.5), voxel_width_*((double)y + 0.5), h-l)); // get heightfield relative to ground
    }
  }
  std::cout << "undercroft removed = " << count << " out of " << heightfield_.rows()*heightfield_.cols() << std::endl;
  drawHeightField("highfield.png", heightfield_);
  drawHeightField("lowfield.png", lowfield_);

  const double curvature_height = 8.0;
#define PARABOLOID
#if defined PARABOLOID
  const double max_diameter_per_height = 1.2; 
  const double min_diameter_per_height = 0.15;

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
  if (verbose)
  {
    std::cout << "num points " << mesh.vertices().size() << std::endl;
    mesh.reduce();
    std::cout << "num points2 " << mesh.vertices().size() << std::endl;
    Mesh mesh2 = mesh;
    for (auto &point: mesh2.vertices())
    {
      double ground_height = lowfield_(int(point[0]/voxel_width_), int(point[1]/voxel_width_));
      point += Eigen::Vector3d(min_bounds_[0], min_bounds_[1], ground_height);
    }
    // ideally we now colour the vertices based on which cluster they're in....

    writePlyMesh("peak_points_mesh.ply", mesh2);
  }
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
  std::vector< std::vector<int> > point_clusters;
  agglomerate(mesh.vertices(), mesh.indexList(), min_diameter_per_height, max_diameter_per_height, point_clusters);
  std::cout << "number found: " << point_clusters.size() << std::endl;
  std::vector<Eigen::Vector3d> &verts = mesh.vertices();

  if (verbose)
  {
    std::vector<Eigen::Vector3d> cloud_points;
    std::vector<double> times;
    std::vector<RGBA> colours;
    for (auto &cluster: point_clusters)
    {
      RGBA colour;
      colour.red = uint8_t(rand()%255);
      colour.green = uint8_t(rand()%255);
      colour.blue = uint8_t(rand()%255);
      colour.alpha = 255;
      for (auto &i: cluster)
      {
        double ground_height = lowfield_(int(verts[i][0]/voxel_width_), int(verts[i][1]/voxel_width));
        Eigen::Vector3d pos = verts[i];
        pos += Eigen::Vector3d(min_bounds_[0], min_bounds_[1], ground_height + 0.05);
        cloud_points.push_back(pos);
        times.push_back(0.0);
        colours.push_back(colour);
      }
    }
    writePlyPointCloud("clusters.ply", cloud_points, times, colours);
  }

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
    tree.tree_tip = min_bounds_ + tip + Eigen::Vector3d(0,0,tree.ground_height);
    tree.radius = tip[2] * max_diameter_per_height * 0.5;
    tree.curvature = -curvature_height/tip[2];
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
    Eigen::Vector3d base = result.tree_tip;
    base[2] = result.ground_height;
    const double tree_radius_to_trunk_radius = 1.0/20.0; // TODO: temporary until we have a better parameter choice
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << result.radius*tree_radius_to_trunk_radius << std::endl;
  }
  return true;
}
}
