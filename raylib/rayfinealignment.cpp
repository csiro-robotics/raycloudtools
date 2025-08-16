// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayfinealignment.h"
#include <nabo/nabo.h>

namespace ray
{

// Convert the set of points into a covariance matrix, and from that into surfel information, using an
// eigendecomposition
void getSurfel(const std::vector<Eigen::Vector3d> &points, const std::vector<int> &ids, Eigen::Vector3d &centroid,
               Eigen::Vector3d &width, Eigen::Matrix3d &mat)
{
  Eigen::Vector3d total(0, 0, 0);
  for (auto &id : ids) total += points[id];
  centroid = total / (double)ids.size();

  Eigen::Matrix3d scatter;
  scatter.setZero();
  for (auto &id : ids)
  {
    Eigen::Vector3d offset = points[id] - centroid;
    scatter += offset * offset.transpose();
  }
  scatter / (double)ids.size();

  // eigendecomposition:
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
  ASSERT(eigen_solver.info() == Eigen::ComputationInfo::Success);
  width = ray::maxVector(eigen_solver.eigenvalues(), Eigen::Vector3d(1e-5, 1e-5, 1e-5));
  // ellipsoid radii are the square root because it is the decomposition of a covariance matrix
  width = Eigen::Vector3d(sqrt(width[0]), sqrt(width[1]), sqrt(width[2]));
  mat = eigen_solver.eigenvectors();
  if (mat.determinant() < 0.0)
    mat.col(0) = -mat.col(0);  // make right-handed, so that we can convert to a quaternion for rendering
}

// Convert clouds_[] into sets of surfels.
void FineAlignment::generateSurfels()
{
  double avg_max_spacing = 0.0;
  double min_spacing = 0.0, max_spacing = 0.0;
  const double min_spacing_scale = 2.0;
  const double max_spacing_scale = 10.0;
  for (int c = 0; c < 2; c++)
  {
    double point_spacing = clouds_[c].estimatePointSpacing();
    ASSERT(point_spacing >= 0.0);
    min_spacing += 0.5 * min_spacing_scale * point_spacing;
    max_spacing += 0.5 * max_spacing_scale * point_spacing;
    avg_max_spacing += 0.5 * max_spacing;
  }
  if (verbose_)
    std::cout << "fine alignment min voxel size: " << min_spacing << "m and maximum voxel size: " << max_spacing
              << "m" << std::endl;

  for (int c = 0; c < 2; c++)
  {
    // 1. decimate quite fine
    std::vector<int64_t> decimated;
    ray::voxelSubsample(clouds_[c].ends, min_spacing, decimated);
    std::vector<Eigen::Vector3d> decimated_points;
    decimated_points.reserve(decimated.size());
    std::vector<Eigen::Vector3d> decimated_starts;
    decimated_starts.reserve(decimated.size());
    centres_[c].setZero();
    for (size_t i = 0; i < decimated.size(); i++)
    {
      if (clouds_[c].rayBounded((int)decimated[i]))
      {
        decimated_points.push_back(clouds_[c].ends[decimated[i]]);
        centres_[c] += decimated_points.back();
        decimated_starts.push_back(clouds_[c].starts[decimated[i]]);
      }
    }
    centres_[c] /= (double)decimated_points.size();

    // 2. find the coarser random candidate points. We just want a fairly even spread but not the voxel centres
    std::vector<int64_t> candidates;
    ray::voxelSubsample(decimated_points, max_spacing, candidates);
    std::vector<Eigen::Vector3d> candidate_points(candidates.size());
    std::vector<Eigen::Vector3d> candidate_starts(candidates.size());
    for (int64_t i = 0; i < (int64_t)candidates.size(); i++)
    {
      candidate_points[i] = decimated_points[candidates[i]];
      candidate_starts[i] = decimated_starts[candidates[i]];
    }

    // Now find all the finely decimated points that are close neighbours of each coarse candidate point
    size_t q_size = candidates.size();
    size_t p_size = decimated_points.size();
    const int search_size = std::min(20, (int)p_size - 1);
    Nabo::NNSearchD *nns;
    Eigen::MatrixXd points_q(3, q_size);
    for (size_t i = 0; i < q_size; i++) 
      points_q.col(i) = candidate_points[i];
    Eigen::MatrixXd points_p(3, p_size);
    for (size_t i = 0; i < p_size; i++) 
      points_p.col(i) = decimated_points[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

    // Run the search
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    nns->knn(points_q, indices, dists2, search_size, 0.01 * max_spacing, 0, max_spacing);
    delete nns;

    // Convert these set of nearest neighbours into surfels
    surfels_[c].reserve(q_size);
    std::vector<int> ids;
    ids.reserve(search_size);
    const size_t min_points_per_ellipsoid = 5;
    for (size_t i = 0; i < q_size; i++)
    {
      ids.clear();
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++) ids.push_back(indices(j, i));
      if (ids.size() < min_points_per_ellipsoid)  // not dense enough
        continue;

      Eigen::Vector3d centroid;
      Eigen::Vector3d width;
      Eigen::Matrix3d mat;
      getSurfel(decimated_points, ids, centroid, width, mat);
      double q1 = width[0] / width[1];
      double q2 = width[1] / width[2];
      if (q2 < q1)  // cylindrical
      {
        if (q2 > 0.5)  // not cylinderical enough
          continue;
        // register two ellipsoids as the normal is ambiguous
        surfels_[c].push_back(Surfel(centroid, mat, width, mat.col(2), false));
        if (c == 1)
          surfels_[c].push_back(Surfel(centroid, mat, width, -mat.col(2), false));
      }
      else  // planar
      {
        Eigen::Vector3d normal = mat.col(0);
        if ((centroid - candidate_starts[i]).dot(normal) > 0.0)
          normal = -normal;
        // now repeat but removing back facing points. This deals better with double walls, which are quite common
        if (!no_normals_)
        {
          for (int j = (int)ids.size() - 1; j >= 0; j--)
          {
            int id = ids[j];
            if ((decimated_points[id] - decimated_starts[id]).dot(normal) > 0.0)
            {
              ids[j] = ids.back();
              ids.pop_back();
            }
          }
        }
        if (ids.size() < min_points_per_ellipsoid)  // not dense enough
          continue;
        getSurfel(decimated_points, ids, centroid, width, mat);
        normal = mat.col(0);
        double q1 = width[0] / width[1];

        if (q1 > 0.5)  // not planar enough
          continue;
        if ((centroid - candidate_starts[i]).dot(normal) > 0.0)
          normal = -normal;
        surfels_[c].push_back(Surfel(centroid, mat, width, normal, true));
        if (c == 1)
          surfels_[c].push_back(Surfel(centroid, mat, width, -normal, true));
      }
    }
  }
  translation_weight_ = 0.4 / avg_max_spacing;  // smaller finds matches further away
}

// Match surfels_[0] to surfels_[1] based on proximity, normal difference and whether it is a plane or cylinder
void FineAlignment::generateSurfelMatches(std::vector<Match> &matches)
{
  std::vector<Eigen::Vector3d> line_starts;
  std::vector<Eigen::Vector3d> line_ends;
  int search_size = 1;
  size_t q_size = surfels_[0].size();
  size_t p_size = surfels_[1].size();
  Nabo::NNSearchD *nns;
  Eigen::MatrixXd points_q(7, q_size);
  for (size_t i = 0; i < q_size; i++)
  {
    Surfel &s = surfels_[0][i];
    Eigen::Vector3d p = s.centroid * translation_weight_;
    p[2] *= 2.0;  // doen't make much difference...
    points_q.col(i) << p, s.normal, s.is_plane ? 1.0 : 0.0;
  }
  Eigen::MatrixXd points_p(7, p_size);
  for (size_t i = 0; i < p_size; i++)
  {
    Surfel &s = surfels_[1][i];
    Eigen::Vector3d p = s.centroid * translation_weight_;
    p[2] *= 2.0;
    points_p.col(i) << p, s.normal, s.is_plane ? 1.0 : 0.0;
  }
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 7);

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, q_size);
  dists2.resize(search_size, q_size);
  nns->knn(points_q, indices, dists2, search_size, ray::kNearestNeighbourEpsilon * max_normal_difference_, 0,
           max_normal_difference_);
  delete nns;

  for (int i = 0; i < (int)q_size; i++)
  {
    for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++)
    {
      Match match;
      match.ids[0] = i;
      match.ids[1] = indices(j, i);
      Surfel &s0 = surfels_[0][i];
      Surfel &s1 = surfels_[1][indices(j, i)];
      if (s0.is_plane != s1.is_plane)
        continue;
      Eigen::Vector3d mid_norm = (s0.normal + s1.normal).normalized();
      if (s0.is_plane)
      {
        match.normal = mid_norm;
        matches.push_back(match);
      }
      else
      {
        // a cylinder is like two normal constraints
        match.normal = mid_norm.cross(Eigen::Vector3d(1, 2, 3)).normalized();
        matches.push_back(match);
        match.normal = mid_norm.cross(match.normal);
        matches.push_back(match);
      }
      line_starts.push_back(s0.centroid);
      line_ends.push_back(s1.centroid);
    }
  }
}

// Convert the correspondences into a linear system to solve
void FineAlignment::buildLinearSystem(const std::vector<Match> &matches, double d, FineAlignment::LinearSystem &system)
{
  // don't go above 30*... or below 10*...
  double square_error = 0.0;
  for (size_t i = 0; i < matches.size(); i++)
  {
    auto &match = matches[i];
    Surfel &s0 = surfels_[0][match.ids[0]];
    Surfel &s1 = surfels_[1][match.ids[1]];
    Eigen::Vector3d positions[2] = { s0.centroid, s1.centroid };
    double error = (positions[1] - positions[0]).dot(match.normal);  // mahabolonis instead?
    double error_sqr;
    if (s0.is_plane)
      error_sqr = ray::sqr(error * translation_weight_);
    else
    {
      Eigen::Vector3d flat = positions[1] - positions[0];
      Eigen::Vector3d norm = s0.normal;
      flat -= norm * flat.dot(norm);
      error_sqr = (flat * translation_weight_).squaredNorm();
    }
    // the normal difference is part of the error,
    if (!no_normals_)
      error_sqr += (s0.normal - s1.normal).squaredNorm();
    double weight = pow(std::max(1.0 - error_sqr / ray::sqr(max_normal_difference_), 0.0), d * d);
    square_error += ray::sqr(error);
    Eigen::Matrix<double, 1, LinearSystem::state_size> a;  // the Jacobian
    a.setZero();

    for (int i = 0; i < 3; i++)  // change in error with change in raycloud translation
      a[i] = match.normal[i];
    for (int i = 0; i < 3; i++)  // change in error with change in raycloud orientation
    {
      Eigen::Vector3d axis(0, 0, 0);
      axis[i] = 1.0;
      a[3 + i] = -(positions[0].cross(axis)).dot(match.normal);
    }
    if (non_rigid_)
    {
      positions[0] -= centres_[0];
      positions[1] -= centres_[1];
      a[6] = ray::sqr(positions[0][0]) * match.normal[0];
      a[7] = ray::sqr(positions[0][0]) * match.normal[1];
      a[8] = ray::sqr(positions[0][1]) * match.normal[0];
      a[9] = ray::sqr(positions[0][1]) * match.normal[1];
      a[10] = positions[0][0] * positions[0][1] * match.normal[0];
      a[11] = positions[0][0] * positions[0][1] * match.normal[1];
    }
    system.At_A += a.transpose() * weight * a;
    system.At_b += a.transpose() * weight * error;
  }
  if (verbose_)
    std::cout << "rmse: " << sqrt(square_error / (double)matches.size()) << std::endl;
}

Eigen::Matrix<double, FineAlignment::LinearSystem::state_size, 1> FineAlignment::LinearSystem::solve(bool verbose)
{
  Eigen::Matrix<double, FineAlignment::LinearSystem::state_size, 1> x = At_A.ldlt().solve(At_b);
  if (verbose)
    std::cout << "least squares shift: " << x[0] << ", " << x[1] << ", " << x[2] << ", rotation: " << x[3] << ", "
              << x[4] << ", " << x[5] << std::endl;
  return x;
}

// Update clouds_[0] and surfels_[0] from the specified transformation
void FineAlignment::updateLinearSystem(std::vector<Match> &matches, const QuadraticTransformation &trans)
{
  Pose shift = trans.getEuclideanPart();
  for (size_t i = 0; i < surfels_[0].size(); i++)
  {
    Eigen::Vector3d &pos = surfels_[0][i].centroid;
    Eigen::Vector3d relPos = pos - centres_[0];
    if (non_rigid_)
      pos += trans.a * ray::sqr(relPos[0]) + trans.b * ray::sqr(relPos[1]) + trans.c * relPos[0] * relPos[1];
    pos = shift * pos;
    surfels_[0][i].normal = shift.rotation * surfels_[0][i].normal;
  }
  Eigen::Quaterniond half_rot(Eigen::AngleAxisd(trans.rotation.norm() / 2.0, trans.rotation.normalized()));
  for (auto &match : matches) 
    match.normal = half_rot * match.normal;

  // NOTE: transforming the whole cloud each time is a bit slow,
  // we should be able to concatenate these transforms and only apply them once at the end
  for (auto &end : clouds_[0].ends)
  {
    Eigen::Vector3d relPos = end - centres_[0];
    if (non_rigid_)
      end += trans.a * ray::sqr(relPos[0]) + trans.b * ray::sqr(relPos[1]) + trans.c * relPos[0] * relPos[1];
    end = shift * end;
  }
}

// The fine grained alignment method
void FineAlignment::align()
{
  // missing rays identified by start and end points being equal
  int n0 = (int)clouds_[0].ends.size()/2;
  int n1 = (int)clouds_[1].ends.size()/2;
  if (clouds_[0].ends[n0] == clouds_[0].starts[n0] || clouds_[1].ends[n1] == clouds_[1].starts[n1])
  {
    std::cout << "Warning: at least some end and start points are coincident, so aligning using two-sided surfaces." << std::endl;
    no_normals_ = true;
  }
  // For each cloud: decimate the cloud to make it even, but still quite detailed, e.g. one point per cubic 10cm
  // Decimate again to pick one point per cubic 1m (for instance)
  // Now match the closest X points in 1 to those in 2, and generate surfel per point in 2.
  generateSurfels();

  // Iteratively reweighted least squares. Iteration loop:
  int max_iterations = 8;
  for (int it = 0; it < max_iterations; it++)
  {
    // Match surfels in cloud0 to those in cloud1
    std::vector<Match> matches;
    generateSurfelMatches(matches);

    // Convert the match constraints into a linear system
    LinearSystem system;
    double d = 20.0 * (double)it / (double)max_iterations;
    buildLinearSystem(matches, d, system);

    // Solve as a weighted least squares problem, to find the transformation of best fit
    QuadraticTransformation perturbation(system.solve(verbose_));

    // Update the ray cloud and surfels from on the transformation of best fit
    updateLinearSystem(matches, perturbation);
  }
}

}  // namespace ray
