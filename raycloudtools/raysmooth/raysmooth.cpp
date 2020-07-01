// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Smooth a ray cloud. Nearby off-surface points are moved onto the nearest surface." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysmooth raycloud" << std::endl;
  exit(exit_code);
}

void smoothPointCloud(std::vector<Eigen::Vector3d> &positions, std::vector<Eigen::Vector3d> &normals, int num_neighbors,
                      int smoothing_iterations, double r_bar)
{
  ASSERT(positions.size() == normals.size());
  ASSERT(num_neighbors > 0);
  ASSERT(num_neighbors <= (int)positions.size());

  std::cout << "smooth_pointcloud with " << positions.size() << " points, " << num_neighbors << " neighbours, "
       << smoothing_iterations << " iters, rbar " << r_bar << std::endl;

  // Set up structures for search (pNumDims,numPoints)
  Eigen::MatrixXd data(6, positions.size());
  for (unsigned int i = 0; i < positions.size(); i++) data.col(i) << positions[i], normals[i];

  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", std::min<size_t>(positions.size(), 8));
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(data, num_neighbors, 0, params);

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(num_neighbors, positions.size());
  dists2.resize(num_neighbors, positions.size());
  nns->knn(data, indices, dists2, num_neighbors, ray::kNearestNeighbourEpsilon, 0);

  // Set up data structures for output

  std::vector<Eigen::Vector3d> smooth_normals = normals;

  const double rbar2 = ray::sqr(r_bar);

  for (int iter = 1; iter < smoothing_iterations; ++iter)
  {
#pragma omp parallel for schedule(guided, 128)
    for (unsigned int i = 0; i < positions.size(); ++i)
    {
      Eigen::Vector3d normal = normals[i];
      Eigen::Matrix3d scatter = normal * normal.transpose();

      for (int j = 0; j < num_neighbors && indices(j, i) > -1; ++j)
      {
        int k = indices(j, i);
        double d = 1.0 - (normals[k].dot(normal));
        double weight = (d > 1.0) ? 0.0 : (1.0 / (1 + ray::sqr(d) / rbar2));
        scatter += weight * normals[k] * normals[k].transpose();
      }

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
      ASSERT(eigen_solver.info() == Success);
      Eigen::Matrix3d eigen_vector = eigen_solver.eigenvectors();

      for (int j = 0; j < 3; ++j) smooth_normals[i][j] = eigen_vector.coeff(j, 2);

      // make sure normal doesn't flip
      if (normal.dot(eigen_vector.col(2)) < 0)
        smooth_normals[i] = -smooth_normals[i];
    }

    // copy over the updated normals
    normals = smooth_normals;
  }

  const double surface_r_bar = 0.05;
  std::vector<Eigen::Vector3d> smooth_points(positions.size());

#pragma omp parallel for schedule(guided, 128)
  for (unsigned int i = 0; i < positions.size(); ++i)
  {
    int j, k;
    Eigen::Vector3d normal = normals[i];
    double t, t0;
    t = t0 = normal.dot(positions[i]);

    for (int iter = 0; iter < 3; ++iter)
    {
      double total_distance = 0;
      double total_weight = 1.0;
      for (j = 0; j < num_neighbors && indices(j, i) > -1; ++j)
      {
        k = indices(j, i);
        if (normal.dot(normals[k]) < cos(45.0 / 180.0 * ray::kPi))
          continue;
        double distance = normal.dot(positions[k]) - t;
        double weight = 1.0 / (1.0 + ray::sqr(distance / surface_r_bar));
        total_distance += weight * distance;
        total_weight += weight;
      }
      t += total_distance / total_weight;
    }
    smooth_points[i] = positions[i] + normal * (t - t0);
  }

  delete nns;
  positions = smooth_points;
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 2)
    usage();

  std::string file = argv[1];
  ray::Cloud cloud;
  cloud.load(file);

  std::vector<Eigen::Vector3d> normals = cloud.generateNormals();

  smoothPointCloud(cloud.ends, normals, 15, 10, 10);

  cloud.save(file.substr(0, file.length() - 4) + "_smooth.ply");

  return true;
}
