// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayclusters.h"
#include <nabo/nabo.h>
#include "raylib/raydebugdraw.h"

namespace ray
{
void clustersAgglomerate(const std::vector<Eigen::Vector3d> &points, double min_diameter, double max_diameter,
                         std::vector<std::vector<int>> &point_clusters)
{
  // 1. get nearest neighbours for each point
  const int search_size = std::min(8, (int)points.size() - 1);
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++) points_p.col(i) = points[i];
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, min_diameter);

  struct Nd
  {
    Nd(int id1, int id2, double dist2)
      : id1(id1)
      , id2(id2)
      , dist2(dist2)
    {}
    int id1, id2;
    double dist2;
  };
  std::vector<Nd> nds;
  for (size_t i = 0; i < points.size(); i++)
  {
    for (int j = 0; j < search_size && indices(j, i) > -1; j++)
    {
      nds.push_back(Nd((int)i, indices(j, i), dists2(j, i)));
    }
  }
  std::sort(nds.begin(), nds.end(), [](const Nd &nd1, const Nd &nd2) { return nd1.dist2 < nd2.dist2; });
  struct Cluster
  {
    Eigen::Vector3d min_bound, max_bound;
    std::vector<int> ids;
    bool active;
  };
  std::vector<Cluster> clusters(points.size());
  std::vector<int> cluster_ids(points.size());
  std::vector<bool> visited(points.size());
  for (size_t i = 0; i < points.size(); i++)
  {
    clusters[i].min_bound = clusters[i].max_bound = points[i];
    clusters[i].active = true;
    clusters[i].ids.push_back((int)i);
    cluster_ids[i] = (int)i;
    visited[i] = false;
  }

  // 2. for each node in turn, from smallest to highest distance, agglomerate
  for (auto &node : nds)
  {
    if (cluster_ids[node.id1] == cluster_ids[node.id2])  // already part of same cluster
      continue;
    int cl1 = cluster_ids[node.id1];
    int cl2 = cluster_ids[node.id2];
    Eigen::Vector3d minb = minVector(clusters[cl1].min_bound, clusters[cl2].min_bound);
    Eigen::Vector3d maxb = maxVector(clusters[cl1].max_bound, clusters[cl2].max_bound);
    Eigen::Vector3d dims = maxb - minb;
    double diam = std::max(dims[0], std::max(dims[1], dims[2]));
    if (diam < max_diameter)  // then merge
    {
      int first = std::min(cl1, cl2);
      int last = std::max(cl1, cl2);
      clusters[first].min_bound = minb;
      clusters[first].max_bound = maxb;
      clusters[first].ids.insert(clusters[first].ids.begin(), clusters[last].ids.begin(), clusters[last].ids.end());
      for (auto &id : clusters[last].ids) cluster_ids[id] = first;
      clusters[last].active = false;
    }
  }
  for (auto &cluster : clusters)
  {
    if (cluster.active)
    {
      point_clusters.push_back(cluster.ids);
    }
  }
}

void clustersKMeans(const std::vector<Eigen::Vector3d> &points, double max_diameter,
                    std::vector<std::vector<int>> &point_clusters)
{
  const double big = 1e10;
  int max_k = 4;
  double best_cluster_diam = big;
  for (int k = 1; k <= max_k; k++)
  {
    // try 2 random starts for robustness
    const int max_trials = 2 * k - 1;
    for (int trial = 1; trial <= max_trials; trial++)
    {
      std::vector<Eigen::Vector3d> seeds(k);
      std::vector<int> ids(points.size());
      for (int i = 0; i < (int)points.size(); i++) ids[i] = i;
      for (int s = 0; s < k; s++)
      {
        int seed = std::rand() % (int)ids.size();
        ids[seed] = ids.back();
        ids.pop_back();
        seeds[s] = points[seed];
      }
      // now do iterative clustering on these k seeds:
      const int max_iterations = 5;
      std::vector<std::vector<int>> clusters(k);
      bool bad = false;
      for (int it = 0; it < max_iterations; it++)
      {
        std::vector<Eigen::Vector3d> means(k, Eigen::Vector3d(0, 0, 0));
        std::vector<double> weights(k, 0.0);
        for (size_t pid = 0; pid < points.size(); pid++)
        {
          const Eigen::Vector3d &point = points[pid];
          double min_dist = big;
          int min_i = -1;
          for (int i = 0; i < k; i++)
          {
            double dist = (point - seeds[i]).squaredNorm();
            if (dist < min_dist)
            {
              min_dist = dist;
              min_i = i;
            }
          }
          means[min_i] += point;
          weights[min_i]++;
          if (it == max_iterations - 1)
            clusters[min_i].push_back((int)pid);
        }
        for (int i = 0; i < k; i++)
        {
          if (weights[i] == 0)
          {
            bad = true;
            break;
          }
          seeds[i] = means[i] / weights[i];
        }
        if (bad)
        {
          std::cout << "bad: somehow a mean is farther than any other mean, skipping to next trial" << std::endl;
          break;
        }
      }
      if (bad)
        continue;
      // now collect the points for each cluster and find their diameter
      double max_diam = -1.0;
      for (int i = 0; i < k; i++)
      {
        Eigen::Vector3d min_bound(big, big, big), max_bound(-big, -big, -big);
        for (auto &id : clusters[i])
        {
          min_bound = minVector(min_bound, points[id]);
          max_bound = maxVector(max_bound, points[id]);
        }
        Eigen::Vector3d dims = max_bound - min_bound;
        double diam = std::max(dims[0], std::max(dims[1], dims[2]));
        max_diam = std::max(max_diam, diam);
      }
      if (max_diam < best_cluster_diam)
      {
        best_cluster_diam = max_diam;
        point_clusters = clusters;
      }
      if (max_diam < max_diameter)
      {
        k = max_k + 1;
        trial = max_trials + 1;  // end the loops
      }
    }
  }
}

std::vector<std::vector<int>> generateClusters(const std::vector<Eigen::Vector3d> &points, double min_diameter,
                                               double max_diameter, bool agglomerate, bool verbose)
{
  std::vector<std::vector<int>> point_clusters;
  // corner cases
  if (points.size() == 1)
    point_clusters.push_back(std::vector<int>(1, 0));
  if (points.size() <= 1)
    return point_clusters;

  if (agglomerate)
    clustersAgglomerate(points, min_diameter, max_diameter, point_clusters);
  else
    clustersKMeans(points, max_diameter, point_clusters);

  if (verbose)
  {
    // now render the points by cluster id
    std::vector<double> shades;
    std::vector<Eigen::Vector3d> ps;
    for (size_t i = 0; i < point_clusters.size(); i++)
    {
      double shade = (double)i / (double)(point_clusters.size() - 1);
      for (auto &id : point_clusters[i])
      {
        shades.push_back(shade);
        ps.push_back(points[id]);
      }
    }
    DebugDraw::instance()->drawCloud(ps, shades, 0);
  }

  return point_clusters;
}


}  // namespace ray
