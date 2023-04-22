// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"

#include <nabo/nabo.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Remove noise from ray clouds. In particular edge noise and isolated point noise." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydenoise raycloud 4 cm     - removes rays that contact more than 4 cm from any other," << std::endl;
  std::cout << "raydenoise raycloud 3 sigmas - removes points more than 3 sigmas from nearest points" << std::endl;
  std::cout << "                    range 4 cm - remove mixed-signal noise that occurs at a range gap." << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayDenoise(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::DoubleArgument sigmas(0.0, 100.0);
  ray::DoubleArgument vox_width(1.0, 100.0);
  ray::TextArgument range_text("range");
  ray::DoubleArgument range(1.0, 1000.0);
  ray::TextArgument cm_text("cm");
  ray::ValueKeyChoice quantity({ &vox_width, &sigmas, &range }, { "cm", "sigmas" });

  bool standard_format = ray::parseCommandLine(argc, argv, { &cloud_file, &quantity });
  bool range_noise = ray::parseCommandLine(argc, argv, { &cloud_file, &range_text, &range, &cm_text });
  if (!standard_format && !range_noise)
    usage();

  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name()))
    usage();

  ray::Cloud new_cloud;
  if (range_noise)  // range-based distance measure. For mixed-points where lidar has contacted two surfaces.
  {
    double range_distance = 0.01 * range.value();
    new_cloud.starts.reserve(cloud.starts.size());
    new_cloud.ends.reserve(cloud.ends.size());
    new_cloud.times.reserve(cloud.times.size());
    new_cloud.colours.reserve(cloud.colours.size());
    // Firstly look at adjacent rays by range. We don't want to throw away large changes,
    // instead, the intermediate of 3 adjacent ranges that is too far from both ends...
    for (int i = 1; i < (int)cloud.starts.size() - 1; i++)
    {
      double range0 = (cloud.ends[i - 1] - cloud.starts[i - 1]).norm();
      double range1 = (cloud.ends[i] - cloud.starts[i]).norm();
      double range2 = (cloud.ends[i + 1] - cloud.starts[i + 1]).norm();
      double min_dist =
        std::min(std::abs(range0 - range2), std::min(std::abs(range1 - range0), std::abs(range2 - range1)));
      if (!cloud.rayBounded(i) || min_dist < range_distance)
        new_cloud.addRay(cloud, i);
    }
    std::cout << cloud.starts.size() - new_cloud.starts.size() << " rays removed with range gaps > "
              << range_distance * 100.0 << " cm." << std::endl;
  }
  else if (quantity.selectedKey() == "cm")  // absolute distance measure
  {
    double distance = 0.01 * vox_width.value();
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;

    Nabo::NNSearchD *nns;
    Nabo::Parameters params("bucketSize", 8);
    std::vector<Eigen::Vector3d> &points = cloud.ends;
    Eigen::MatrixXd points_p(3, points.size());
    for (unsigned int i = 0; i < points.size(); i++) 
    {
      points_p.col(i) = points[i];
    }
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

    // Run the search
    const int search_size = 1;
    indices.resize(search_size, points.size());
    dists2.resize(search_size, points.size());
    nns->knn(points_p, indices, dists2, search_size, ray::kNearestNeighbourEpsilon, 0);
    indices.resize(0, 0);
    delete nns;

    new_cloud.starts.reserve(cloud.starts.size());
    new_cloud.ends.reserve(cloud.ends.size());
    new_cloud.times.reserve(cloud.times.size());
    new_cloud.colours.reserve(cloud.colours.size());
    for (int i = 0; i < (int)points.size(); i++)
    {
      if (!cloud.rayBounded(i) || (dists2(0, i) < 1e10 && dists2(0, i) < ray::sqr(distance)))
        new_cloud.addRay(cloud, i);
    }
    std::cout << cloud.starts.size() - new_cloud.starts.size() << " rays removed with ends further than "
              << distance * 100.0 << " cm from any other." << std::endl;
  }
  else if (quantity.selectedKey() == "sigmas")  // scale-invariant distance measure. Same as Mahalanobis distance
  {
    std::vector<Eigen::Vector3d> centroids;
    std::vector<Eigen::Vector3d> dimensions;
    std::vector<Eigen::Matrix3d> matrices;
    Eigen::MatrixXi indices;

    const int search_size = std::min(10, (int)cloud.ends.size() - 1);
    cloud.getSurfels(search_size, &centroids, nullptr, &dimensions, &matrices, &indices);

    new_cloud.starts.reserve(cloud.starts.size());
    new_cloud.ends.reserve(cloud.ends.size());
    new_cloud.times.reserve(cloud.times.size());
    new_cloud.colours.reserve(cloud.colours.size());
    Eigen::Vector3d dims(0, 0, 0);
    double cnt = 0.0;
    double nums = 0;
    for (size_t i = 0; i < matrices.size(); i++)
    {
      bool is_noise = false;
      if (cloud.rayBounded(i))
      {
        if (indices(0, i) == Nabo::NNSearchD::InvalidIndex)  // no neighbours in range, we consider this as noise
          continue;
        int other_i = indices(0, i);
        Eigen::Vector3d vec = cloud.ends[i] - centroids[other_i];
        Eigen::Vector3d newVec = matrices[other_i].transpose() * vec;
        newVec[0] /= dimensions[other_i][0];
        newVec[1] /= dimensions[other_i][1];
        newVec[2] /= dimensions[other_i][2];
        int num = 0;
        for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++) num = j + 1;
        nums += (double)num;
        dims += dimensions[other_i];
        cnt++;
        double scale2 = newVec.squaredNorm();
        is_noise = scale2 > sigmas.value() * sigmas.value();
      }
      if (!is_noise)
        new_cloud.addRay(cloud, i);
    }
    dims /= cnt;
    std::cout << "average dimensions: " << dims.transpose() << ", average num neighbours: " << nums / cnt << std::endl;
    std::cout << cloud.starts.size() - new_cloud.starts.size()
              << " rays removed with nearest neighbour sigma more than " << sigmas.value() << std::endl;
  }

  new_cloud.save(cloud_file.nameStub() + "_denoised.ply");
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDenoise, argc, argv);
}