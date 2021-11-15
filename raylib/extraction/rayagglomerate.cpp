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
static const double wrap_gradient = 1.0;

void Forest::agglomerate(const std::vector<Eigen::Vector3d> &points, const std::vector<Eigen::Vector3i> index_list, double min_diameter_per_height, double max_diameter_per_height, std::vector<Cluster> &point_clusters)
{
  struct Nd
  {
    Nd(int id1, int id2, double dist) : id1(id1), id2(id2), dist(dist) {}
    int id1, id2;
    double dist;
  };
  std::vector<Nd> nds;
  for (auto &ind: index_list) // doubles up on edges, which slows down the sort, but is safely discarded in later code
  {
    for (int i = 0; i<3; i++)
    {
      int id1 = ind[i];
      int id2 = ind[(i+1)%3];
      Eigen::Vector3d diff = points[id1]-points[id2];
      double height = std::abs(diff[2]);
      diff[2] = 0.0;

      double dist = diff.norm();
      double h = height / (wrap_gradient * dist);
      if (h > 1.01)
        std::cout << "how did this happen?" << std::endl;
      double height_penalty_factor = 0.5; // 0 to 1. Make apparent distance less when there is a height difference, because large height differences we cannot see as much downward slope
      dist *= 1.0 - h*height_penalty_factor;

      nds.push_back(Nd(std::min(id1, id2), std::max(id1, id2), dist));
    }
  }
  std::sort(nds.begin(), nds.end(), [](const Nd &nd1, const Nd &nd2){ return nd1.dist < nd2.dist; });

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
  int num_no_peaks = 0;
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
    const double search_up_gradient = 0.2;
    if (closest_i >= 0 && std::sqrt(min_dist2) < search_up_gradient*points[closest_i][2])
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
      num_no_peaks++;
      // now how to we relay this information??
      Eigen::Vector3d tip = trunk.first - min_bounds_;
      tip[2] = maxheight;
      Cluster new_cluster;
      new_cluster.min_bound = new_cluster.max_bound = tip;
      new_cluster.trunk_id = c;
      clusters.push_back(new_cluster);
    }
  }
  std::cout << num_no_peaks << " out of " << trunks_.size() << " trunks have no mesh peaks above them, canopy found directly above trunks instead" << std::endl;

  // 2. for each node in turn, from smallest to highest distance, agglomerate
  for (auto &node: nds)
  {
    int cl1 = cluster_ids[node.id1];
    int cl2 = cluster_ids[node.id2];
    if (cl1 == cl2) // already part of same cluster
      continue; 
    if (clusters[cl1].trunk_id >=0 && clusters[cl2].trunk_id >=0 && clusters[cl1].trunk_id != clusters[cl2].trunk_id) // don't merge from different trunks
    {
  //    std::cout << "clusters " << cl1 << " and " << cl2 << "cannot merge between two different trunks " << clusters[cl1].trunk_id << " and " << clusters[cl2].trunk_id << std::endl;
      continue;
    } 
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = 1.0 + std::max(dims[0], dims[1]);
    double growth_ratio = 1.0;

    double mean_height = (minb[2] + maxb[2])/2.0; 

    double mid_node_height = (points[node.id1][2]+points[node.id2][2])/2.0; // TODO: divide by this mid height before sorting
    if (node.dist > mid_node_height * min_diameter_per_height / growth_ratio) // /growth_ratio means you separate trees more easily when connectinng two distinct blobs
    {
      // about to keep as isolated cluster. So check if there is space underneath, and only isolate if there is
      #define MERGE_IF_NO_SPACE // haven't seen this make a difference yet on viewed data... but it ought to be a worthwhile condition
      #if defined MERGE_IF_NO_SPACE 
      Eigen::Vector3d tip;
      bool space1 = findSpace(clusters[cl1], points, tip);
      bool space2 = findSpace(clusters[cl2], points, tip);
      if (space1 && space2) // space is found for each so they're separate, the end
      #endif
        continue; 
    }
    if (diam < mean_height * max_diameter_per_height / growth_ratio) // then merge
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

void Forest::renderAgglomeration(const std::vector<Cluster> &point_clusters, const std::vector<Eigen::Vector3d> &verts, const std::string &cloud_name_stub)
{
  std::vector<Eigen::Vector3d> cloud_points;
  std::vector<double> times;
  std::vector<RGBA> colours;
  RGBA colour;
  colour.alpha = 255;
  for (auto &cluster: point_clusters)
  {
    if (cluster.ids.empty())
      continue;
    colour.red = uint8_t(rand()%255);
    colour.green = uint8_t(rand()%255);
    colour.blue = uint8_t(rand()%255);
    for (auto &i: cluster.ids)
    {
      double ground_height = lowfield_(int(verts[i][0]/voxel_width_), int(verts[i][1]/voxel_width_));
      Eigen::Vector3d pos = verts[i];
      pos += Eigen::Vector3d(min_bounds_[0], min_bounds_[1], ground_height + 0.05);
      cloud_points.push_back(pos);
      times.push_back(0.0);
      colours.push_back(colour);
    }
    Eigen::Vector3d tip;
    bool found = findSpace(cluster, verts, tip);
    double height_per_radius = 50.0;
    double rad = tip[2] / height_per_radius;
    int x = int(tip[0]/voxel_width_);
    int y = int(tip[1]/voxel_width_);
    if (x < 0 || x >= lowfield_.rows() || y < 0 || y >= lowfield_.cols())
    {
      std::cout << "bad lookup: " << x << ", " << y << std::endl;
      continue;
    }
    tip[2] = lowfield_(x, y);
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
        if (cloud_points.size() == 29798)
        {
          std::cout << "third: " << Eigen::Vector3d(x, y, height).transpose() << std::endl;
        }        
        cloud_points.push_back(Eigen::Vector3d(x, y, height));
        if (cloud_points.back().norm() > 1e5)
          std::cout << " bad space: " << cloud_points.back() << std::endl;
        times.push_back(0.0);
        colour.red = colour.green = colour.blue = (uint8_t)(255.0*spacefield_(i,j));
        colours.push_back(colour);
      }
    }
  }

  writePlyPointCloud(cloud_name_stub + "_clusters.ply", cloud_points, times, colours);
}

} // namespace ray