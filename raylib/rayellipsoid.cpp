// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tom Lowe, Kazys Stepanas
#include "rayellipsoid.h"

#include "raycloud.h"
#include "rayprogress.h"

#include <nabo/nabo.h>

#include <memory>

#if RAYLIB_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB

namespace ray
{
void generateEllipsoids(std::vector<Ellipsoid> *ellipsoids, Eigen::Vector3d *bounds_min, Eigen::Vector3d *bounds_max,
                        const Cloud &cloud, Progress *progress)
{
  ellipsoids->clear();
  ellipsoids->resize(cloud.rayCount());
  const int search_size = std::min(16, (int)cloud.rayCount() - 1);
  const double max_double = std::numeric_limits<double>::max();
  Eigen::Vector3d ellipsoids_min(max_double, max_double, max_double);
  Eigen::Vector3d ellipsoids_max(-max_double, -max_double, -max_double);
  Nabo::Parameters params("bucketSize", 8);

  if (progress)
  {
    progress->begin("generateEllipsoids - KDTree", 2);
  }

  Eigen::MatrixXd points_p(3, cloud.ends.size());
  for (size_t i = 0; i < cloud.rayCount(); ++i)
  {
    points_p.col(i) = cloud.ends[i];
  }
  std::unique_ptr<Nabo::NNSearchD> nns(Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3));

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, cloud.rayCount());
  dists2.resize(search_size, cloud.rayCount());

  if (progress)
  {
    progress->increment();
  }
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
  nns.reset(nullptr);

  if (progress)
  {
    progress->increment();
    progress->end();
    progress->begin("generateEllipsoids", cloud.ends.size());
  }
  const auto generate_ellipsoid = [&](size_t i)  //
  {
    Ellipsoid &ellipsoid = (*ellipsoids)[i];
    ellipsoid.clear();
    ellipsoid.transient = false;
    ellipsoid.opacity = 1.0;

    // Increment progress here as we have multiple exit points
    if (progress)
    {
      progress->increment();
    }

    if (!cloud.rayBounded(i))
    {
      return;
    }

    Eigen::Matrix3d scatter;
    scatter.setZero();
    Eigen::Vector3d centroid(0, 0, 0);
    double num_neighbours = 0;
    for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; ++j)
    {
      int index = indices(j, i);
      if (cloud.rayBounded(index))
      {
        centroid += cloud.ends[index];
        num_neighbours++;
      }
    }
    if (num_neighbours < 4)
    {
      return;
    }
    centroid /= num_neighbours;
    for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
    {
      int index = indices(j, i);
      if (cloud.rayBounded(index))
      {
        Eigen::Vector3d offset = cloud.ends[index] - centroid;
        scatter += offset * offset.transpose();
      }
    }
    scatter /= num_neighbours;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Eigen::ComputationInfo::Success);

    Eigen::Vector3d eigen_value = eigen_solver.eigenvalues();
    Eigen::Matrix3d eigen_vector = eigen_solver.eigenvectors();

    ellipsoid.pos = centroid;
    double scale = 1.7;  // this scale roughly matches the dimensions of a uniformly dense ellipsoid
    eigen_value[0] = scale * sqrt(std::max(1e-10, eigen_value[0]));
    eigen_value[1] = scale * sqrt(std::max(1e-10, eigen_value[1]));
    eigen_value[2] = scale * sqrt(std::max(1e-10, eigen_value[2]));
    ellipsoid.eigen_mat.row(0) = (eigen_vector.col(0) / eigen_value[0]).cast<float>();
    ellipsoid.eigen_mat.row(1) = (eigen_vector.col(1) / eigen_value[1]).cast<float>();
    ellipsoid.eigen_mat.row(2) = (eigen_vector.col(2) / eigen_value[2]).cast<float>();
    ellipsoid.time = cloud.times[i];
    ellipsoid.setExtents(eigen_vector, eigen_value);
  };

#if RAYLIB_WITH_TBB
  tbb::parallel_for<size_t>(0, cloud.rayCount(), generate_ellipsoid);

  for (size_t i = 0; i < ellipsoids->size(); ++i)
  {
    Ellipsoid &ellipsoid = (*ellipsoids)[i];
    const auto ellipsoid_min = ellipsoid.pos - ellipsoid.extents.cast<double>();
    const auto ellipsoid_max = ellipsoid.pos + ellipsoid.extents.cast<double>();

    ellipsoids_min.x() = std::min(ellipsoids_min.x(), ellipsoid_min.x());
    ellipsoids_min.y() = std::min(ellipsoids_min.y(), ellipsoid_min.y());
    ellipsoids_min.z() = std::min(ellipsoids_min.z(), ellipsoid_min.z());
    ellipsoids_max.x() = std::max(ellipsoids_max.x(), ellipsoid_max.x());
    ellipsoids_max.y() = std::max(ellipsoids_max.y(), ellipsoid_max.y());
    ellipsoids_max.z() = std::max(ellipsoids_max.z(), ellipsoid_max.z());
  }
#else   // RAYLIB_WITH_TBB
  const size_t count = cloud.rayCount();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < count; ++i)
  {
    generate_ellipsoid(i);
  }
  for (size_t i = 0; i < ellipsoids->size(); ++i)
  {
    Ellipsoid &ellipsoid = (*ellipsoids)[i];
    const auto ellipsoid_min = ellipsoid.pos - ellipsoid.extents.cast<double>();
    const auto ellipsoid_max = ellipsoid.pos + ellipsoid.extents.cast<double>();

    ellipsoids_min.x() = std::min(ellipsoids_min.x(), ellipsoid_min.x());
    ellipsoids_min.y() = std::min(ellipsoids_min.y(), ellipsoid_min.y());
    ellipsoids_min.z() = std::min(ellipsoids_min.z(), ellipsoid_min.z());
    ellipsoids_max.x() = std::max(ellipsoids_max.x(), ellipsoid_max.x());
    ellipsoids_max.y() = std::max(ellipsoids_max.y(), ellipsoid_max.y());
    ellipsoids_max.z() = std::max(ellipsoids_max.z(), ellipsoid_max.z());
  }
#endif  // RAYLIB_WITH_TBB

  if (bounds_min)
  {
    *bounds_min = ellipsoids_min;
  }

  if (bounds_max)
  {
    *bounds_max = ellipsoids_max;
  }
}
}  // namespace ray
