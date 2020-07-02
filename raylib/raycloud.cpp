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

#include <iostream>
#include <limits>
#include <set>

namespace ray
{

void Cloud::clear()
{
  starts.clear();
  ends.clear();
  times.clear();
  colours.clear();
}

void Cloud::save(const std::string &file_name) const
{
  std::string name = file_name;
  if (name.substr(name.length() - 4) != ".ply")
    name += ".ply";
  writePly(name, starts, ends, times, colours);
}

bool Cloud::load(const std::string &file_name)
{
  // look first for the raycloud PLY
  if (file_name.substr(file_name.size() - 4) == ".ply")
    return loadPLY(file_name);
  if (std::ifstream((file_name + ".ply").c_str(), std::ios::in))
    return loadPLY(file_name + ".ply");

  // otherwise, look for a .laz and _traj.txt file by that name
  if (std::ifstream((file_name + ".laz").c_str(), std::ios::in) && 
      std::ifstream((file_name + "_traj.txt").c_str(), std::ios::in))
    return loadLazTraj(file_name + ".laz", file_name + "_traj.txt");

  return false;
}

bool Cloud::load(const std::string &point_cloud, const std::string &traj_file)
{
  std::string name_end = point_cloud.substr(point_cloud.size() - 4);
  if (name_end == ".ply")
    readPly(point_cloud, starts, ends, times, colours, false); // special case of reading a non-ray-cloud ply
  else if (name_end == ".laz" || name_end == ".las")
    readLas(point_cloud, ends, times, colours, 1);
  else
  {
    std::cout << "Error converting unknown type: " << point_cloud << std::endl;
    return false;
  }

  Trajectory trajectory;
  trajectory.load(traj_file);
  calculateStarts(trajectory);
  return true;
}

bool Cloud::loadPLY(const std::string &file)
{
  return readPly(file, starts, ends, times, colours, true);
}

bool Cloud::loadLazTraj(const std::string &laz_file, const std::string &traj_file)
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
    std::cout << "can only recalculate when a trajectory is available" << std::endl;
}

Eigen::Vector3d Cloud::calcMinBound() const
{
  Eigen::Vector3d min_v(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      min_v = minVector(min_v, minVector(starts[i], ends[i]));
  }
  return min_v;
}

Eigen::Vector3d Cloud::calcMaxBound() const
{
  Eigen::Vector3d max_v(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
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
    progress->begin("calcBounds", rayCount());
  }

  *min_bounds = Eigen::Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
  *max_bounds = Eigen::Vector3d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
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
  std::vector<int> valids;
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

void Cloud::decimate(double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> *voxel_set)
{
  std::vector<int64_t> subsample = voxelSubsample(ends, voxel_width, voxel_set);
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

void Cloud::getSurfels(int search_size, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                       std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats, 
                       Eigen::MatrixXi *neighbour_indices)
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
  std::vector<int> ray_ids;
  ray_ids.reserve(ends.size());
  for (unsigned int i = 0; i < ends.size(); i++)
    if (rayBounded(i))
      ray_ids.push_back(i);
  Eigen::MatrixXd points_p(3, ray_ids.size());
  for (unsigned int i = 0; i < ray_ids.size(); i++) 
    points_p.col(i) = ends[ray_ids[i]];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, ray_ids.size());
  dists2.resize(search_size, ray_ids.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
  delete nns;

  if (neighbour_indices)
    neighbour_indices->resize(search_size, ends.size());
  for (int i = 0; i < (int)ray_ids.size(); i++)
  {
    int ray_id = ray_ids[i];
    if (neighbour_indices)
    {
      int j;
      for (j = 0; j < search_size && indices(j, i) > -1; j++) 
        (*neighbour_indices)(j, ray_id) = ray_ids[indices(j, i)];
      if (j < search_size)
        (*neighbour_indices)(j, ray_id) = -1;
    }

    Eigen::Vector3d centroid = ends[ray_id];
    int num;
    for (num = 0; num < search_size && indices(num, i) > -1; num++) centroid += ends[ray_ids[indices(num, i)]];
    centroid /= (double)(num + 1);
    if (centroids)
      (*centroids)[ray_id] = centroid;

    Eigen::Matrix3d scatter = (ends[ray_id] - centroid) * (ends[ray_id] - centroid).transpose();
    for (int j = 0; j < num; j++)
    {
      Eigen::Vector3d offset = ends[ray_ids[indices(j, i)]] - centroid;
      scatter += offset * offset.transpose();
    }
    scatter /= (double)(num + 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Success);
    if (normals)
    {
      Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
      if ((ends[ray_id] - starts[ray_id]).dot(normal) > 0.0)
        normal = -normal;
      (*normals)[ray_id] = normal;
    }
    if (dimensions)
    {
      Eigen::Vector3d eigenvals = maxVector(Eigen::Vector3d(1e-10, 1e-10, 1e-10), eigen_solver.eigenvalues());
      (*dimensions)[ray_id] = Eigen::Vector3d(std::sqrt(eigenvals[0]), std::sqrt(eigenvals[1]), std::sqrt(eigenvals[2]));
    }
    if (mats)
      (*mats)[ray_id] = eigen_solver.eigenvectors();
  }
}

// starts are required to get the normal the right way around
std::vector<Eigen::Vector3d> Cloud::generateNormals(int search_size)
{
  std::vector<Eigen::Vector3d> normals;
  getSurfels(search_size, NULL, &normals, NULL, NULL, NULL);
  return normals;
}

double Cloud::estimatePointSpacing() const
{
  // two-iteration estimation, modelling the point distribution by the below exponent.
  // larger exponents (towards 2.5) match thick forests, lower exponents (towards 2) match smooth terrain and surfaces
  const double cloud_exponent = 2.0; // model num_points = (cloud_width/voxel_width)^cloud_exponent

  Eigen::Vector3d min_bound, max_bound;
  calcBounds(&min_bound, &max_bound, kBFEnd);
  Eigen::Vector3d extent = max_bound - min_bound;
  int num_points = 0;
  for (unsigned int i = 0; i < ends.size(); i++)
    if (rayBounded(i))
      num_points++;
  double cloud_width = pow(extent[0]*extent[1]*extent[2], 1.0/3.0); // an average
  double voxel_width = cloud_width / pow((double)ends.size(), 1.0/cloud_exponent);
  voxel_width *= 5.0; // we want to use a larger width because this process only works when the width is an overestimation
  std::cout << "initial voxel width estimate: " << voxel_width << std::endl;
  double num_voxels = 0;
  std::set<Eigen::Vector3i, Vector3iLess> test_set;
  for (unsigned int i = 0; i < ends.size(); i++)
  {
    if (rayBounded(i))
    {
      const Eigen::Vector3d &point = ends[i];
      Eigen::Vector3i place(int(std::floor(point[0] / voxel_width)), int(std::floor(point[1] / voxel_width)),
                            int(std::floor(point[2] / voxel_width)));
      if (test_set.find(place) == test_set.end())
      {
        test_set.insert(place);
        num_voxels++;
      }
    }
  }
  double points_per_voxel = (double)num_points / num_voxels;
  double width = voxel_width / pow(points_per_voxel, 1.0/cloud_exponent);
  std::cout << "estimated point spacing: " << width << std::endl;
  return width;
}

void Cloud::split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr)
{
  for (int i = 0; i < (int)ends.size(); i++)
  {
    Cloud &cloud = fptr(i) ? cloud2 : cloud1;
    cloud.addRay(*this, i);
  }
}

void Cloud::addRay(const Eigen::Vector3d &start, const Eigen::Vector3d &end, double time, const RGBA &colour)
{
  starts.push_back(start);
  ends.push_back(end);
  times.push_back(time);
  colours.push_back(colour);
}


void Cloud::addRay(const Cloud &other_cloud, size_t index)
{
  starts.push_back(other_cloud.starts[index]);
  ends.push_back(other_cloud.ends[index]);
  times.push_back(other_cloud.times[index]);
  colours.push_back(other_cloud.colours[index]);
}

} // namespace ray
