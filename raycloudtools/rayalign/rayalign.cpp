// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raylib/raycloud.h"
#include "raylib/rayalignment.h"
#include "raylib/raydebugdraw.h"

#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <complex>

void usage(int exit_code = 0)
{
  std::cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayalign raycloudA raycloudB" << std::endl;
  std::cout << "                             --nonrigid    - nonrigid (quadratic) alignment" << std::endl;
  std::cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << std::endl;
  std::cout
    << "                             --local    - fine alignment only, assumes clouds are already approximately aligned"
    << std::endl;
  exit(exit_code);
}

void getSurfel(const std::vector<Eigen::Vector3d> &points, const std::vector<int> &ids, Eigen::Vector3d &centroid, Eigen::Vector3d &width,
               Eigen::Matrix3d &mat)
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

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
  ASSERT(eigen_solver.info() == Success);
  width = ray::maxVector(eigen_solver.eigenvalues(), Eigen::Vector3d(1e-5, 1e-5, 1e-5));
  width = Eigen::Vector3d(sqrt(width[0]), sqrt(width[1]), sqrt(width[2]));
  mat = eigen_solver.eigenvectors();
}

int main(int argc, char *argv[])
{
  // TODO: This method works when there is more than 30% overlap.
  // For less, we can additionally repeat this procedure for small scale cubic sections (perhaps 20% of map width)
  // then use the fft power spectrum (low res) as a descriptor into the knn to find matching locations,
  // then pick the best match.
  ray::DebugDraw &draw = *ray::DebugDraw::init(argc, argv, "rayalign");

  if (argc < 3 || argc > 6)
    usage();
  bool verbose = false;
  bool non_rigid = false;
  bool local_only = false;
  for (int a = 3; a < argc; a++)
  {
    if (std::string(argv[a]) == "--verbose" || std::string(argv[a]) == "-v")
      verbose = true;
    else if (std::string(argv[a]) == "--nonrigid" || std::string(argv[a]) == "-n")
      non_rigid = true;
    else if (std::string(argv[a]) == "--local" || std::string(argv[a]) == "-l")
      local_only = true;
    else
      usage();
  }

  std::string file_a = argv[1];
  std::string file_b = argv[2];
  std::string file_stub = file_a;
  if (file_stub.substr(file_stub.length() - 4) == ".ply")
    file_stub = file_stub.substr(0, file_stub.length() - 4);

  ray::AlignTranslationYaw aligner;
  aligner.clouds[0].load(file_a);
  aligner.clouds[1].load(file_b);

  if (!local_only)
  {
    aligner.alignCloud0ToCloud1(0.5, verbose);
    if (verbose)
      aligner.clouds[0].save(file_stub + "_coarse_aligned.ply");
  }

  // Next the fine grained alignment.
  // Method:
  // 1. for each cloud: decimate the cloud to make it even, but still quite detailed, e.g. one point per cubic 10cm
  // 2. decimate again to pick one point per cubic 1m (for instance)
  // 3. now match the closest X points in 1 to those in 2, and generate surfel per point in 2.
  // 4. repeatedly: match surfels in cloud0 to those in cloud1
  // 5. solve a weighted least squares to find the transform of best fit
  // 6. apply the transform to cloud0 and save it out.

  std::vector<Eigen::Matrix3d> matrices[2];

  std::vector<Eigen::Vector3d> centroids[2];
  std::vector<Eigen::Vector3d> normals[2];
  std::vector<Eigen::Vector3d> widths[2];
  std::vector<bool> is_plane[2];
  Eigen::Vector3d centres[2];
  double point_spacings[2];
  double avg_max_spacing = 0.0; 
  for (int c = 0; c < 2; c++)
  {
    point_spacings[c] = aligner.clouds[c].estimatePointSpacing();
    const double min_spacing_scale = 2.0;
    const double max_spacing_scale = 20.0;
    double min_spacing = min_spacing_scale*point_spacings[c];
    double max_spacing = max_spacing_scale*point_spacings[c];
    avg_max_spacing += 0.5*max_spacing;
    std::cout << "fine alignment min voxel size: " << min_spacing << "m and maximum voxel size: " << max_spacing << "m" << std::endl;
    // 1. decimate quite fine
    std::vector<int64_t> decimated = ray::voxelSubsample(aligner.clouds[c].ends, min_spacing);
    std::vector<Eigen::Vector3d> decimated_points;
    decimated_points.reserve(decimated.size());
    std::vector<Eigen::Vector3d> decimated_starts;
    decimated_starts.reserve(decimated.size());
    centres[c].setZero();
    for (size_t i = 0; i < decimated.size(); i++)
    {
      if (aligner.clouds[c].rayBounded((int)decimated[i]))
      {
        decimated_points.push_back(aligner.clouds[c].ends[decimated[i]]);
        centres[c] += decimated_points.back();
        decimated_starts.push_back(aligner.clouds[c].starts[decimated[i]]);
      }
    }
    centres[c] /= (double)decimated_points.size(); 

    // 2. find the coarser random candidate points. We just want a fairly even spread but not the voxel centres
    std::vector<int64_t> candidates = ray::voxelSubsample(decimated_points, max_spacing);
    std::vector<Eigen::Vector3d> candidate_points(candidates.size());
    std::vector<Eigen::Vector3d> candidate_starts(candidates.size());
    for (int64_t i = 0; i < (int64_t)candidates.size(); i++)
    {
      candidate_points[i] = decimated_points[candidates[i]];
      candidate_starts[i] = decimated_starts[candidates[i]];
    }

    int search_size = 20;
    size_t q_size = candidates.size();
    size_t p_size = decimated_points.size();

    Nabo::NNSearchD *nns;
    Eigen::MatrixXd points_q(3, q_size);
    for (size_t i = 0; i < q_size; i++) points_q.col(i) = candidate_points[i];
    Eigen::MatrixXd points_p(3, p_size);
    for (size_t i = 0; i < p_size; i++) points_p.col(i) = decimated_points[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

    // Run the search
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    nns->knn(points_q, indices, dists2, search_size, 0.01*max_spacing, 0, max_spacing);
    delete nns;

    centroids[c].reserve(q_size);
    normals[c].reserve(q_size);
    widths[c].reserve(q_size);
    matrices[c].reserve(q_size);
    is_plane[c].reserve(q_size);
    std::vector<int> ids;
    ids.reserve(search_size);
    const size_t min_points_per_ellipsoid = 5;
    for (size_t i = 0; i < q_size; i++)
    {
      ids.clear();
      for (int j = 0; j < search_size && indices(j, i) > -1; j++) 
        ids.push_back(indices(j, i));
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
        matrices[c].push_back(mat);
        centroids[c].push_back(centroid);
        normals[c].push_back(mat.col(2));
        widths[c].push_back(width);
        is_plane[c].push_back(false);
        if (c == 1)
        {
          matrices[c].push_back(mat);
          centroids[c].push_back(centroid);
          normals[c].push_back(-mat.col(2));
          widths[c].push_back(width);
          is_plane[c].push_back(false);
        }
      }
      else
      {
        Eigen::Vector3d normal = mat.col(0);
        if ((centroid - candidate_starts[i]).dot(normal) > 0.0)
          normal = -normal;
        // now repeat but removing back facing points. This deals better with double walls, which are quite common
        for (int j = (int)ids.size() - 1; j >= 0; j--)
        {
          int id = ids[j];
          if ((decimated_points[id] - decimated_starts[id]).dot(normal) > 0.0)
          {
            ids[j] = ids.back();
            ids.pop_back();
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
        matrices[c].push_back(mat);
        centroids[c].push_back(centroid);
        normals[c].push_back(normal);
        widths[c].push_back(width);
        is_plane[c].push_back(true);
      }
    }
    draw.drawCloud(decimated_points, 0.5 + 0.4 * (double)c, c);
  }
  draw.drawEllipsoids(centroids[1], matrices[1], widths[1], Eigen::Vector3d(0, 1, 0), 1);

  // Now we have the ellipsoids, we need to do a nearest neighbour on this set of data, trying to match orientation as
  // well as location
  struct Match
  {
    int ids[2];
    Eigen::Vector3d normal;
  };
  std::vector<Match> matches;
  {
    // with these matches we now run the iterative reweighted least squares..
    ray::Pose pose = ray::Pose::identity();
    int max_iterations = 8;
    const double translation_weight = 0.4 / avg_max_spacing;  // smaller finds matches further away
    const double max_normal_difference = 0.5;
    for (int it = 0; it < max_iterations; it++)
    {
      std::vector<Eigen::Vector3d> line_starts;
      std::vector<Eigen::Vector3d> line_ends;
      int search_size = 1;
      size_t q_size = centroids[0].size();
      size_t p_size = centroids[1].size();
      Nabo::NNSearchD *nns;
      Eigen::MatrixXd points_q(7, q_size);
      for (size_t i = 0; i < q_size; i++)
      {
        Eigen::Vector3d p = centroids[0][i] * translation_weight;
        p[2] *= 2.0;  // doen't make much difference...
        points_q.col(i) << p, normals[0][i], is_plane[0][i] ? 1.0 : 0.0;
      }
      Eigen::MatrixXd points_p(7, p_size);
      for (size_t i = 0; i < p_size; i++)
      {
        Eigen::Vector3d p = centroids[1][i] * translation_weight;
        p[2] *= 2.0;
        points_p.col(i) << p, normals[1][i], is_plane[1][i] ? 1.0 : 0.0;
      }
      nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 7);

      // Run the search
      Eigen::MatrixXi indices;
      Eigen::MatrixXd dists2;
      indices.resize(search_size, q_size);
      dists2.resize(search_size, q_size);
      nns->knn(points_q, indices, dists2, search_size, 0.01*max_normal_difference, 0, max_normal_difference);
      delete nns;
      matches.clear();

      for (int i = 0; i < (int)q_size; i++)
      {
        for (int j = 0; j < search_size && indices(j, i) > -1; j++)
        {
          Match match;
          match.ids[0] = i;
          match.ids[1] = indices(j, i);
          bool plane = is_plane[0][i];
          if (plane != is_plane[1][indices(j, i)])
            continue;
          Eigen::Vector3d mid_norm = (normals[0][i] + normals[1][indices(j, i)]).normalized();
          if (plane)
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
          line_starts.push_back(centroids[0][i]);
          line_ends.push_back(centroids[1][indices(j, i)]);
        }
      }
      draw.drawLines(line_starts, line_ends);
      draw.drawEllipsoids(centroids[0], matrices[0], widths[0], Eigen::Vector3d(1, 0, 0), 0);

      // don't go above 30*... or below 10*...
      double d = 20.0 * (double)it / (double)max_iterations;
      const int state_size = 12;
      Eigen::Matrix<double, state_size, state_size> at_a;
      at_a.setZero();
      Eigen::Matrix<double, state_size, 1> at_b;
      at_b.setZero();
      double square_error = 0.0;
      for (size_t i = 0; i < matches.size(); i++)
      {
        auto &match = matches[i];
        Eigen::Vector3d positions[2] = { centroids[0][match.ids[0]], centroids[1][match.ids[1]] };
        double error = (positions[1] - positions[0]).dot(match.normal);  // mahabolonis instead?
        double error_sqr;
        if (is_plane[0][match.ids[0]])
          error_sqr = ray::sqr(error * translation_weight);
        else
        {
          Eigen::Vector3d flat = positions[1] - positions[0];
          Eigen::Vector3d norm = normals[0][match.ids[0]];
          flat -= norm * flat.dot(norm);
          error_sqr = (flat * translation_weight).squaredNorm();
        }
        // the normal difference is part of the error,
        error_sqr += (normals[0][match.ids[0]] - normals[1][match.ids[1]]).squaredNorm();
        double weight = pow(std::max(1.0 - error_sqr / ray::sqr(max_normal_difference), 0.0), d * d);
        square_error += ray::sqr(error);
        Eigen::Matrix<double, state_size, 1> at;  // the Jacobian
        at.setZero();

        for (int i = 0; i < 3; i++)  // change in error with change in raycloud translation
          at[i] = match.normal[i];
        for (int i = 0; i < 3; i++)  // change in error with change in raycloud orientation
        {
          Eigen::Vector3d axis(0, 0, 0);
          axis[i] = 1.0;
          at[3 + i] = -(positions[0].cross(axis)).dot(match.normal);
        }
        if (non_rigid) 
        {
          positions[0] -= centres[0];
          positions[1] -= centres[1];
          at[6] = ray::sqr(positions[0][0]) * match.normal[0];
          at[7] = ray::sqr(positions[0][0]) * match.normal[1];
          at[8] = ray::sqr(positions[0][1]) * match.normal[0];
          at[9] = ray::sqr(positions[0][1]) * match.normal[1];
          at[10] = positions[0][0] * positions[0][1] * match.normal[0];
          at[11] = positions[0][0] * positions[0][1] * match.normal[1];
        }
        at_a += at * weight * at.transpose();
        at_b += at * weight * error;
      }
      Eigen::Matrix<double, state_size, 1> x = at_a.ldlt().solve(at_b);
      std::cout << "rmse: " << sqrt(square_error / (double)matches.size()) << std::endl;
      std::cout << "least squares shift: " << x[0] << ", " << x[1] << ", " << x[2] << ", rotation: " << x[3] << ", " << x[4]
           << ", " << x[5] << std::endl;
      Eigen::Vector3d rot(x[3], x[4], x[5]);
      Eigen::Vector3d a(x[6], x[7], 0), b(x[8], x[9], 0), c(x[10], x[11], 0);
      double angle = rot.norm();
      rot.normalize();
      ray::Pose shift(Eigen::Vector3d(x[0], x[1], x[2]), Eigen::Quaterniond(Eigen::AngleAxisd(angle, rot)));
      pose = shift * pose;
      for (size_t i = 0; i < centroids[0].size(); i++)
      {
        Eigen::Vector3d &pos = centroids[0][i];
        Eigen::Vector3d relPos = pos - centres[0];
        if (non_rigid)
          pos += a * ray::sqr(relPos[0]) + b * ray::sqr(relPos[1]) + c * relPos[0] * relPos[1];
        pos = shift * pos;
        normals[0][i] = shift.rotation * normals[0][i];
      }
      Eigen::Quaterniond half_rot(Eigen::AngleAxisd(angle / 2.0, rot));
      for (auto &match : matches) match.normal = half_rot * match.normal;

      // TODO: transforming the whole cloud each time is a bit slow,
      // we should be able to concatenate these transforms and only apply them once at the end
      for (auto &end : aligner.clouds[0].ends)
      {
        Eigen::Vector3d relPos = end - centres[0];
        if (non_rigid)
          end += a * ray::sqr(relPos[0]) + b * ray::sqr(relPos[1]) + c * relPos[0] * relPos[1];
        end = shift * end;
      }
    }
  }

  aligner.clouds[0].save(file_stub + "_aligned.ply");
  return true;
}
