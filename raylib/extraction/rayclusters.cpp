// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayclusters.h"
#include <nabo/nabo.h>

namespace ray
{
// take the input points and separate into clusters based on a minimum and maximum separation diameter criterion
// this is a form of agglomerative clustering
void clustersAgglomerate(const std::vector<Eigen::Vector3d> &points, double min_diameter, double max_diameter,
                         std::vector<std::vector<int>> &point_clusters)
{
  // 1. get nearest neighbours for each point
  const int search_size = std::min(8, static_cast<int>(points.size()) - 1);
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++)
  {
    points_p.col(i) = points[i];
  }
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, min_diameter);

  // temporary node structure in order to sort the neighbours by distance
  struct Nd
  {
    Nd(int id1, int id2, double dist2)
      : id1(id1)
      , id2(id2)
      , dist2(dist2)
    {}
    int id1, id2;
    double dist2;  // square distance
  };
  std::vector<Nd> nds;
  for (size_t i = 0; i < points.size(); i++)
  {
    for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
    {
      nds.push_back(Nd(static_cast<int>(i), indices(j, i), dists2(j, i)));
    }
  }
  std::sort(nds.begin(), nds.end(), [](const Nd &nd1, const Nd &nd2) { return nd1.dist2 < nd2.dist2; });

  // temporary cluster structure
  struct Cluster
  {
    Eigen::Vector3d min_bound, max_bound;
    std::vector<int> ids;
    bool active;
  };
  std::vector<Cluster> clusters(points.size());
  std::vector<int> cluster_ids(points.size());
  // initialise the clusters
  for (size_t i = 0; i < points.size(); i++)
  {
    clusters[i].min_bound = clusters[i].max_bound = points[i];
    clusters[i].active = true;
    clusters[i].ids.push_back(static_cast<int>(i));
    cluster_ids[i] = static_cast<int>(i);
  }

  // 2. for each node in turn, from smallest to highest distance, agglomerate
  for (auto &node : nds)
  {
    if (cluster_ids[node.id1] == cluster_ids[node.id2])  // already part of same cluster
    {
      continue;
    }
    const int cl1 = cluster_ids[node.id1];
    const int cl2 = cluster_ids[node.id2];
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = std::max(dims[0], std::max(dims[1], dims[2]));
    if (diam < max_diameter)  // then merge
    {
      const int first = std::min(cl1, cl2);
      const int last = std::max(cl1, cl2);
      clusters[first].min_bound = minb;
      clusters[first].max_bound = maxb;
      clusters[first].ids.insert(clusters[first].ids.begin(), clusters[last].ids.begin(), clusters[last].ids.end());
      for (auto &id : clusters[last].ids)
      {
        cluster_ids[id] = first;
      }
      clusters[last].active = false;
    }
  }
  // convert the vector of clusters into a vector of sets of point indices
  for (auto &cluster : clusters)
  {
    if (cluster.active)
    {
      point_clusters.push_back(cluster.ids);
    }
  }
}

/// generate clusters from the set of points, with optional debug rending of the output
void generateClusters(std::vector<std::vector<int>> &point_clusters, const std::vector<Eigen::Vector3d> &points,
                      double min_diameter, double max_diameter)
{
  // corner cases
  if (points.size() == 1)
  {
    point_clusters.push_back(std::vector<int>(1, 0));
  }
  if (points.size() <= 1)
  {
    return;
  }

  clustersAgglomerate(points, min_diameter, max_diameter, point_clusters);
}


}  // namespace ray
