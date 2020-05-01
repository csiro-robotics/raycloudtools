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

using namespace std;
using namespace Eigen;
using namespace ray;

void usage(int exit_code = 0)
{
  cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << endl;
  cout << "usage:" << endl;
  cout << "rayalign raycloudA raycloudB" << endl;
  cout << "                             --rigid    - rigid alignment only" << endl;
  cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << endl;
  cout
    << "                             --local    - fine alignment only, assumes clouds are already approximately aligned"
    << endl;
  exit(exit_code);
}

void getSurfel(const vector<Vector3d> &points, const vector<int> &ids, Vector3d &centroid, Vector3d &width,
               Matrix3d &mat)
{
  Vector3d total(0, 0, 0);
  for (auto &id : ids) total += points[id];
  centroid = total / (double)ids.size();

  Matrix3d scatter;
  scatter.setZero();
  for (auto &id : ids)
  {
    Vector3d offset = points[id] - centroid;
    scatter += offset * offset.transpose();
  }
  scatter / (double)ids.size();

  SelfAdjointEigenSolver<Matrix3d> eigen_solver(scatter.transpose());
  ASSERT(eigenSolver.info() == Success);
  width = maxVector(eigen_solver.eigenvalues(), Vector3d(1e-5, 1e-5, 1e-5));
  width = Vector3d(sqrt(width[0]), sqrt(width[1]), sqrt(width[2]));
  mat = eigen_solver.eigenvectors();
}

int main(int argc, char *argv[])
{
  // TODO: This method works when there is more than 30% overlap.
  // For less, we can additionally repeat this procedure for small scale cubic sections (perhaps 20% of map width)
  // then use the fft power spectrum (low res) as a descriptor into the knn to find matching locations,
  // then pick the best match.
  DebugDraw &draw = *DebugDraw::init(argc, argv, "rayalign");

  if (argc < 3 || argc > 6)
    usage();
  bool verbose = false;
  bool rigid_only = false;
  bool local_only = false;
  for (int a = 3; a < argc; a++)
  {
    if (string(argv[a]) == "--verbose" || string(argv[a]) == "-v")
      verbose = true;
    else if (string(argv[a]) == "--rigid" || string(argv[a]) == "-r")
      rigid_only = true;
    else if (string(argv[a]) == "--local" || string(argv[a]) == "-l")
      local_only = true;
    else
      usage();
  }

  string file_a = argv[1];
  string file_b = argv[2];
  string file_stub = file_a;
  if (file_stub.substr(file_stub.length() - 4) == ".ply")
    file_stub = file_stub.substr(0, file_stub.length() - 4);

  AlignTranslationYaw aligner;
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

  vector<Matrix3d> matrices[2];

  vector<Vector3d> centroids[2];
  vector<Vector3d> normals[2];
  vector<Vector3d> widths[2];
  vector<bool> is_plane[2];
  for (int c = 0; c < 2; c++)
  {
    // 1. decimate quite fine
    vector<int> decimated = voxelSubsample(aligner.clouds[c].ends, 0.1);
    vector<Vector3d> decimated_points;
    decimated_points.reserve(decimated.size());
    vector<Vector3d> decimated_starts;
    decimated_starts.reserve(decimated.size());
    for (int i = 0; i < (int)decimated.size(); i++)
    {
      if (aligner.clouds[c].rayBounded(decimated[i]))
      {
        decimated_points.push_back(aligner.clouds[c].ends[decimated[i]]);
        decimated_starts.push_back(aligner.clouds[c].starts[decimated[i]]);
      }
    }

    // 2. find the coarser random candidate points. We just want a fairly even spread but not the voxel centres
    vector<int> candidates = voxelSubsample(decimated_points, 1.0);
    vector<Vector3d> candidate_points(candidates.size());
    vector<Vector3d> candidate_starts(candidates.size());
    for (int i = 0; i < (int)candidates.size(); i++)
    {
      candidate_points[i] = decimated_points[candidates[i]];
      candidate_starts[i] = decimated_starts[candidates[i]];
    }

    int search_size = 20;
    size_t q_size = candidates.size();
    size_t p_size = decimated_points.size();

    Nabo::NNSearchD *nns;
    MatrixXd points_q(3, q_size);
    for (size_t i = 0; i < q_size; i++) points_q.col(i) = candidate_points[i];
    MatrixXd points_p(3, p_size);
    for (size_t i = 0; i < p_size; i++) points_p.col(i) = decimated_points[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

    // Run the search
    MatrixXi indices;
    MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    nns->knn(points_q, indices, dists2, search_size, 0.01, 0, 1.0);
    delete nns;

    centroids[c].reserve(q_size);
    normals[c].reserve(q_size);
    widths[c].reserve(q_size);
    matrices[c].reserve(q_size);
    is_plane[c].reserve(q_size);
    vector<int> ids;
    ids.reserve(search_size);
    const size_t min_points_per_ellipsoid = 5;
    for (size_t i = 0; i < q_size; i++)
    {
      ids.clear();
      for (int j = 0; j < search_size && indices(j, i) > -1; j++) ids.push_back(indices(j, i));
      if (ids.size() < min_points_per_ellipsoid)  // not dense enough
        continue;

      Vector3d centroid;
      Vector3d width;
      Matrix3d mat;
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
        Vector3d normal = mat.col(0);
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
  draw.drawEllipsoids(centroids[1], matrices[1], widths[1], Vector3d(0, 1, 0), 1);

  // Now we have the ellipsoids, we need to do a nearest neighbour on this set of data, trying to match orientation as
  // well as location
  struct Match
  {
    int ids[2];
    Vector3d normal;
  };
  vector<Match> matches;
  {
    // with these matches we now run the iterative reweighted least squares..
    Pose pose = Pose::identity();
    int max_iterations = 8;
    const double translation_weight = 0.4;  // smaller finds matches further away
    const double max_normal_difference = 0.5;
    for (int it = 0; it < max_iterations; it++)
    {
      //   if (it < 3) // the distribution of re-matching within the iterations is open to adjustment
      {
        vector<Vector3d> line_starts;
        vector<Vector3d> line_ends;
        int search_size = 1;
        size_t q_size = centroids[0].size();
        size_t p_size = centroids[1].size();
        Nabo::NNSearchD *nns;
        MatrixXd points_q(7, q_size);
        for (size_t i = 0; i < q_size; i++)
        {
          Vector3d p = centroids[0][i] * translation_weight;
          p[2] *= 2.0;  // doen't make much difference...
          points_q.col(i) << p, normals[0][i], is_plane[0][i] ? 1.0 : 0.0;
        }
        MatrixXd points_p(7, p_size);
        for (size_t i = 0; i < p_size; i++)
        {
          Vector3d p = centroids[1][i] * translation_weight;
          p[2] *= 2.0;
          points_p.col(i) << p, normals[1][i], is_plane[1][i] ? 1.0 : 0.0;
        }
        nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 7);

        // Run the search
        MatrixXi indices;
        MatrixXd dists2;
        indices.resize(search_size, q_size);
        dists2.resize(search_size, q_size);
        nns->knn(points_q, indices, dists2, search_size, 0.01, 0, max_normal_difference);
        delete nns;
        matches.clear();

        for (int i = 0; i < (int)q_size; i++)
        {
          // shall we pick only two-way matches? Not for now...
          for (int j = 0; j < search_size && indices(j, i) > -1; j++)
          {
            Match match;
            match.ids[0] = i;
            match.ids[1] = indices(j, i);
            bool plane = is_plane[0][i];
            if (plane != is_plane[1][indices(j, i)])
              continue;
            Vector3d mid_norm = (normals[0][i] + normals[1][indices(j, i)]).normalized();
            if (plane)
            {
              match.normal = mid_norm;
              matches.push_back(match);
            }
            else
            {
              // a cylinder is like two normal constraints
              match.normal = mid_norm.cross(Vector3d(1, 2, 3)).normalized();
              matches.push_back(match);
              match.normal = mid_norm.cross(match.normal);
              matches.push_back(match);
            }
            line_starts.push_back(centroids[0][i]);
            line_ends.push_back(centroids[1][indices(j, i)]);
          }
        }
        draw.drawLines(line_starts, line_ends);
        draw.drawEllipsoids(centroids[0], matrices[0], widths[0], Vector3d(1, 0, 0), 0);
      }

      // don't go above 30*... or below 10*...
      double d = 20.0 * (double)it / (double)max_iterations;
      const int state_size = 12;
      Matrix<double, state_size, state_size> at_a;
      at_a.setZero();
      Matrix<double, state_size, 1> at_b;
      at_b.setZero();
      double square_error = 0.0;
      for (int i = 0; i < (int)matches.size(); i++)
      {
        auto &match = matches[i];
        Vector3d pos[2] = { centroids[0][match.ids[0]], centroids[1][match.ids[1]] };
        double error = (pos[1] - pos[0]).dot(match.normal);  // mahabolonis instead?
        double error_sqr;
        if (is_plane[0][match.ids[0]])
          error_sqr = sqr(error * translation_weight);
        else
        {
          Vector3d flat = pos[1] - pos[0];
          Vector3d norm = normals[0][match.ids[0]];
          flat -= norm * flat.dot(norm);
          error_sqr = (flat * translation_weight).squaredNorm();
        }
        // the normal difference is part of the error,
        error_sqr += (normals[0][match.ids[0]] - normals[1][match.ids[1]]).squaredNorm();
        double weight = pow(max(1.0 - error_sqr / sqr(max_normal_difference), 0.0), d * d);
        square_error += sqr(error);
        Matrix<double, state_size, 1> at;  // the Jacobian
        at.setZero();

        for (int i = 0; i < 3; i++)  // change in error with change in raycloud translation
          at[i] = match.normal[i];
        for (int i = 0; i < 3; i++)  // change in error with change in raycloud orientation
        {
          Vector3d axis(0, 0, 0);
          axis[i] = 1.0;
          at[3 + i] = -(pos[0].cross(axis)).dot(match.normal);
        }
        if (!rigid_only)  // give the aligner a chance to rigidly align first
        {
          at[6] = sqr(pos[0][0]) * match.normal[0];
          at[7] = sqr(pos[0][0]) * match.normal[1];
          at[8] = sqr(pos[0][1]) * match.normal[0];
          at[9] = sqr(pos[0][1]) * match.normal[1];
          at[10] = pos[0][0] * pos[0][1] * match.normal[0];
          at[11] = pos[0][0] * pos[0][1] * match.normal[1];
        }
        at_a += at * weight * at.transpose();
        at_b += at * weight * error;
      }
      Matrix<double, state_size, 1> x = at_a.ldlt().solve(at_b);
      cout << "rmse: " << sqrt(square_error / (double)matches.size()) << endl;
      cout << "least squares shift: " << x[0] << ", " << x[1] << ", " << x[2] << ", rotation: " << x[3] << ", " << x[4]
           << ", " << x[5] << endl;
      Vector3d rot(x[3], x[4], x[5]);
      Vector3d a(x[6], x[7], 0), b(x[8], x[9], 0), c(x[10], x[11], 0);
      double angle = rot.norm();
      rot.normalize();
      Pose shift(Vector3d(x[0], x[1], x[2]), Quaterniond(AngleAxisd(angle, rot)));
      pose = shift * pose;
      for (int i = 0; i < (int)centroids[0].size(); i++)
      {
        Vector3d &pos = centroids[0][i];
        if (!rigid_only)
          pos += a * sqr(pos[0]) + b * sqr(pos[1]) + c * pos[0] * pos[1];
        pos = shift * pos;
        normals[0][i] = shift.rotation * normals[0][i];
      }
      Quaterniond half_rot(AngleAxisd(angle / 2.0, rot));
      for (auto &match : matches) match.normal = half_rot * match.normal;

      // TODO: transforming the whole cloud each time is a bit slow,
      // we should be able to concatenate these transforms and only apply them once at the end
      for (auto &end : aligner.clouds[0].ends)
      {
        if (!rigid_only)
          end += a * sqr(end[0]) + b * sqr(end[1]) + c * end[0] * end[1];
        end = shift * end;
      }
    }
  }

  aligner.clouds[0].save(file_stub + "_aligned.ply");
  return true;
}
