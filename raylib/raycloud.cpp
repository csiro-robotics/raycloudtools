// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"

#include "raydebugdraw.h"
#include "raylaz.h"
#include "rayply.h"
#include "rayprogress.h"
#include "raytrajectory.h"

#include <nabo/nabo.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>

using namespace std;
using namespace Eigen;
using namespace ray;

namespace
{
/// Log a @c std::chrono::clock::duration to an output stream.
///
/// The resulting string displays in the smallest possible unit to show three three
/// decimal places with display units ranging from seconds to nanoseconds. The table below
/// shows some example times.
///
/// Time(s)     | Display
/// ----------- | --------
/// 0.000000018 | 18ns
/// 0.000029318 | 29.318us
/// 0.0295939   | 29.593ms
/// 0.93        | 930ms
/// 15.023      | 15.023s
/// 15.000025   | 15.000s
///
/// Note that times are truncated, not rounded.
///
/// @tparam D The duration type of the form @c std::chrono::clock::duration.
/// @param out The output stream to log to.
/// @param duration The duration to convert to string.
template <typename D>
inline std::ostream &logDuration(std::ostream &out, const D &duration)
{
  const bool negative = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() < 0;
  const char *sign = (!negative) ? "" : "-";
  D abs_duration = (!negative) ? duration : duration * -1;
  auto s = std::chrono::duration_cast<std::chrono::seconds>(abs_duration).count();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(abs_duration).count();
  ms = ms % 1000;

  if (s)
  {
    out << sign << s << "." << std::setw(3) << std::setfill('0') << ms << "s";
  }
  else
  {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(abs_duration).count();
    us = us % 1000;

    if (ms)
    {
      out << sign << ms << "." << std::setw(3) << std::setfill('0') << us << "ms";
    }
    else
    {
      auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(abs_duration).count();
      ns = ns % 1000;

      if (us)
      {
        out << sign << us << "." << std::setw(3) << std::setfill('0') << ns << "us";
      }
      else
      {
        out << sign << ns << "ns";
      }
    }
  }
  return out;
}
}  // namespace

void Cloud::clear()
{
  starts.clear();
  ends.clear();
  times.clear();
  colours.clear();
}

void Cloud::save(const std::string &file_name) const
{
  string name = file_name;
  if (name.substr(name.length() - 4) != ".ply")
    name += ".ply";
  writePly(name, starts, ends, times, colours);
}

bool Cloud::load(const std::string &file_name)
{
  // look first for the raycloud PLY
  if (file_name.substr(file_name.size() - 4) == ".ply")
    return loadPLY(file_name);
  if (ifstream((file_name + ".ply").c_str(), ios::in))
    return loadPLY(file_name + ".ply");

  // otherwise, look for a .laz and _traj.txt file by that name
  if (ifstream((file_name + ".laz").c_str(), ios::in) && ifstream((file_name + "_traj.txt").c_str(), ios::in))
    return loadLazTraj(file_name + ".laz", file_name + "_traj.txt");

  return false;
}

bool Cloud::load(const std::string &point_cloud, const std::string &traj_file)
{
  string name_end = point_cloud.substr(point_cloud.size() - 4);
  if (name_end == ".ply")
    readPly(point_cloud, starts, ends, times, colours);
  else if (name_end == ".laz" || name_end == ".las")
    readLas(point_cloud, ends, times, colours, 1);
  else
  {
    cout << "Error converting unknown type: " << point_cloud << endl;
    return false;
  }

  Trajectory trajectory;
  trajectory.load(traj_file);
  calculateStarts(trajectory);
  return true;
}

bool Cloud::loadPLY(const string &file)
{
  return readPly(file, starts, ends, times, colours);
}

bool Cloud::loadLazTraj(const string &laz_file, const string &traj_file)
{
  bool success = readLas(laz_file, ends, times, colours, 1);
  if (!success)
    return false;
  Trajectory trajectory;
  trajectory.load(traj_file);
  calculateStarts(trajectory);
  return true;
}

void Cloud::calculateStarts(const Trajectory &trajectory)
{
  // Aha!, problem in calculating starts when times are not ordered.
  if (trajectory.nodes.size() > 0)
  {
    int n = 1;
    starts.resize(ends.size());
    for (size_t i = 0; i < ends.size(); i++)
    {
      while ((times[i] > trajectory.nodes[n].time) && n < (int)trajectory.nodes.size() - 1) n++;
      double blend =
        (times[i] - trajectory.nodes[n - 1].time) / (trajectory.nodes[n].time - trajectory.nodes[n - 1].time);
      starts[i] =
        trajectory.nodes[n - 1].pose.position +
        (trajectory.nodes[n].pose.position - trajectory.nodes[n - 1].pose.position) * clamped(blend, 0.0, 1.0);
    }
  }
  else
    cout << "can only recalculate when a trajectory is available" << endl;
}

Vector3d Cloud::calcMinBound() const
{
  Vector3d min_v(numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max());
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      min_v = minVector(min_v, minVector(starts[i], ends[i]));
  }
  return min_v;
}

Vector3d Cloud::calcMaxBound() const
{
  Vector3d max_v(numeric_limits<double>::lowest(), numeric_limits<double>::lowest(), numeric_limits<double>::lowest());
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      max_v = maxVector(max_v, maxVector(starts[i], ends[i]));
  }
  return max_v;
}

bool Cloud::calcBounds(Eigen::Vector3d *min_bounds, Eigen::Vector3d *max_bounds, unsigned flags, Progress *progress) const
{
  if (rayCount() == 0)
  {
    return false;
  }

  if (progress)
  {
    progress->reset("calcBounds", rayCount());
  }

  *min_bounds = Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
  *max_bounds = Vector3d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                         std::numeric_limits<double>::lowest());
  bool invalid_bounds = true;
  for (size_t i = 0; i < rayCount(); ++i)
  {
    if (rayBounded(i))
    {
      invalid_bounds = false;
      if (flags & kBFEnd)
      {
        *min_bounds = minVector(*min_bounds, ends[i]);
        *max_bounds = maxVector(*max_bounds, ends[i]);
      }
      if (flags & kBFStart)
      {
        *min_bounds = minVector(*min_bounds, starts[i]);
        *max_bounds = maxVector(*max_bounds, starts[i]);
      }
    }

    if (progress)
    {
      progress->increment();
    }
  }

  return !invalid_bounds;
}

void Cloud::transform(const Pose &pose, double time_delta)
{
  for (int i = 0; i < (int)starts.size(); i++)
  {
    starts[i] = pose * starts[i];
    ends[i] = pose * ends[i];
    times[i] += time_delta;
  }
}

void Cloud::removeUnboundedRays()
{
  vector<int> valids;
  for (int i = 0; i < (int)ends.size(); i++)
    if (rayBounded(i))
      valids.push_back(i);
  for (int i = 0; i < (int)valids.size(); i++)
  {
    starts[i] = starts[valids[i]];
    ends[i] = ends[valids[i]];
    times[i] = times[valids[i]];
    colours[i] = colours[valids[i]];
  }
  starts.resize(valids.size());
  ends.resize(valids.size());
  times.resize(valids.size());
  colours.resize(valids.size());
}

void Cloud::decimate(double voxel_width)
{
  vector<int64_t> subsample = voxelSubsample(ends, voxel_width);
  for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
  {
    int64_t id = subsample[i];
    starts[i] = starts[id];
    ends[i] = ends[id];
    colours[i] = colours[id];
    times[i] = times[id];
  }
  starts.resize(subsample.size());
  ends.resize(subsample.size());
  colours.resize(subsample.size());
  times.resize(subsample.size());
}

void Cloud::getSurfels(int search_size, vector<Vector3d> *centroids, vector<Vector3d> *normals,
                       vector<Vector3d> *dimensions, vector<Matrix3d> *mats, MatrixXi *neighbour_indices)
{
  // simplest scheme... find 3 nearest neighbours and do cross product
  if (centroids)
    centroids->resize(ends.size());
  if (normals)
    normals->resize(ends.size());
  if (dimensions)
    dimensions->resize(ends.size());
  if (mats)
    mats->resize(ends.size());
  Nabo::NNSearchD *nns;
  vector<int> ray_ids;
  ray_ids.reserve(ends.size());
  for (unsigned int i = 0; i < ends.size(); i++)
    if (rayBounded(i))
      ray_ids.push_back(i);
  MatrixXd points_p(3, ray_ids.size());
  for (unsigned int i = 0; i < ray_ids.size(); i++) points_p.col(i) = ends[ray_ids[i]];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(search_size, ray_ids.size());
  dists2.resize(search_size, ray_ids.size());
  nns->knn(points_p, indices, dists2, search_size, 0.01, 0, 1.0);  // TODO: needs to sort here
  delete nns;

  if (neighbour_indices)
    neighbour_indices->resize(search_size, ends.size());
  for (int i = 0; i < (int)ray_ids.size(); i++)
  {
    int ii = ray_ids[i];
    if (neighbour_indices)
    {
      int j;
      for (j = 0; j < search_size && indices(j, i) > -1; j++) (*neighbour_indices)(j, ii) = ray_ids[indices(j, i)];
      if (j < search_size)
        (*neighbour_indices)(j, ii) = -1;
    }

    Vector3d centroid = ends[i];
    int num;
    for (num = 0; num < search_size && indices(num, i) > -1; num++) centroid += ends[ray_ids[indices(num, i)]];
    centroid /= (double)(num + 1);
    if (centroids)
      (*centroids)[i] = centroid;

    Matrix3d scatter = (ends[i] - centroid) * (ends[i] - centroid).transpose();
    for (int j = 0; j < num; j++)
    {
      Vector3d offset = ends[ray_ids[indices(j, i)]] - centroid;
      scatter += offset * offset.transpose();
    }
    scatter /= (double)(num + 1);

    SelfAdjointEigenSolver<Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Success);
    if (normals)
    {
      Vector3d normal = eigen_solver.eigenvectors().col(0);
      if ((ends[i] - starts[i]).dot(normal) > 0.0)
        normal = -normal;
      (*normals)[i] = normal;
    }
    if (dimensions)
    {
      Vector3d eigenvals = maxVector(Vector3d(1e-10, 1e-10, 1e-10), eigen_solver.eigenvalues());
      (*dimensions)[i] = Vector3d(sqrt(eigenvals[0]), sqrt(eigenvals[1]), sqrt(eigenvals[2]));
    }
    if (mats)
      (*mats)[i] = eigen_solver.eigenvectors();
  }
}

// starts are required to get the normal the right way around
vector<Vector3d> Cloud::generateNormals(int search_size)
{
  vector<Vector3d> normals;
  getSurfels(search_size, NULL, &normals, NULL, NULL, NULL);
  return normals;
}

double Cloud::estimatePointSpacing() const
{
  double v_width = 0.25;
  double num_voxels = 0;
  double num_points = 0;
  std::set<Eigen::Vector3i, Vector3iLess> test_set;
  for (unsigned int i = 0; i < ends.size(); i++)
  {
    if (rayBounded(i))
    {
      num_points++;
      const Vector3d &point = ends[i];
      Eigen::Vector3i place(int(std::floor(point[0] / v_width)), int(std::floor(point[1] / v_width)),
                            int(std::floor(point[2] / v_width)));
      if (test_set.find(place) == test_set.end())
      {
        test_set.insert(place);
        num_voxels++;
      }
    }
  }

  double width =
    0.25 *
    sqrt(num_voxels /
         num_points);  // since points roughly represent 2D surfaces. Also matches empirical tests of optimal speed
  cout << "estimated point spacing: " << width << endl;
  return width;
}

void Cloud::split(Cloud &cloud1, Cloud &cloud2, function<bool(int i)> fptr)
{
  for (int i = 0; i < (int)ends.size(); i++)
  {
    Cloud &cloud = fptr(i) ? cloud2 : cloud1;
    cloud.starts.push_back(starts[i]);
    cloud.ends.push_back(ends[i]);
    cloud.times.push_back(times[i]);
    cloud.colours.push_back(colours[i]);
  }
}
