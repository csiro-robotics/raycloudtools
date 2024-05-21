// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"

#include "raylaz.h"
#include "rayply.h"
#include "rayprogress.h"

#include <nabo/nabo.h>

#include <iostream>
#include <limits>
#include <set>
// #define OUTPUT_CLOUD_MOMENTS // useful for setting up unit tests comparisons

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
  writePlyRayCloud(name, starts, ends, times, colours);
}

bool Cloud::load(const std::string &file_name, bool check_extension, int min_num_rays)
{
  // look first for the raycloud PLY
  if (file_name.substr(file_name.size() - 4) == ".ply" || !check_extension)
    return loadPLY(file_name, min_num_rays);

  std::cerr << "Attempting to load ray cloud " << file_name << " which doesn't have expected file extension .ply"
            << std::endl;
  return false;
}

bool Cloud::loadPLY(const std::string &file, int min_num_rays)
{
  bool res = readPly(file, starts, ends, times, colours, true);
  if ((int)ends.size() < min_num_rays)
    return false;
#if defined OUTPUT_CLOUD_MOMENTS // Only used to supply data to unit tests
  Eigen::Array<double, 22, 1> mom = getMoments();
  std::cout << "stats: " << std::endl;
  for (int i = 0; i < mom.rows(); i++) 
  {
    std::cout << ", " << mom[i];
  }
  std::cout << std::endl;
#endif  // defined OUTPUT_CLOUD_MOMENTS
  return res;
}

Eigen::Vector3d Cloud::calcMinBound() const
{
  Eigen::Vector3d min_v(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max());
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      min_v = minVector(min_v, minVector(starts[i], ends[i]));
  }
  return min_v;
}

Eigen::Vector3d Cloud::calcMaxBound() const
{
  Eigen::Vector3d max_v(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                        std::numeric_limits<double>::lowest());
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      max_v = maxVector(max_v, maxVector(starts[i], ends[i]));
  }
  return max_v;
}

bool Cloud::calcBounds(Eigen::Vector3d *min_bounds, Eigen::Vector3d *max_bounds, unsigned flags,
                       Progress *progress) const
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

void Cloud::decimate(double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> &voxel_set)
{
  std::vector<int64_t> subsample;
  voxelSubsample(ends, voxel_width, subsample, voxel_set);
  for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
  {
    const int64_t id = subsample[i];
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

void Cloud::eigenSolve(const std::vector<int> &ray_ids, const Eigen::MatrixXi &indices, int index, int num_neighbours,
                       Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> &solver, Eigen::Vector3d &centroid) const
{
  int ray_id = ray_ids[index];
  centroid = ends[ray_id];
  for (int j = 0; j < num_neighbours; j++) centroid += ends[ray_ids[indices(j, index)]];
  centroid /= (double)(num_neighbours + 1);
  Eigen::Matrix3d scatter = (ends[ray_id] - centroid) * (ends[ray_id] - centroid).transpose();
  for (int j = 0; j < num_neighbours; j++)
  {
    Eigen::Vector3d offset = ends[ray_ids[indices(j, index)]] - centroid;
    scatter += offset * offset.transpose();
  }
  scatter /= (double)(num_neighbours + 1);
  solver.compute(scatter.transpose());
  ASSERT(solver.info() == Eigen::ComputationInfo::Success);
}

void Cloud::getSurfels(int search_size, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                       std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                       Eigen::MatrixXi *neighbour_indices, double max_distance, bool reject_back_facing_rays) const
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
  for (unsigned int i = 0; i < ray_ids.size(); i++) points_p.col(i) = ends[ray_ids[i]];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, ray_ids.size());
  dists2.resize(search_size, ray_ids.size());
  if (max_distance != 0.0)
    nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, max_distance);
  else
    nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
  delete nns;

  if (neighbour_indices)
  {
    neighbour_indices->resize(search_size, ends.size());
    for (int i = 0; i<neighbour_indices->rows(); i++)
    {
      for (int j = 0; j < neighbour_indices->cols(); j++)
      {
        (*neighbour_indices)(i, j) = -1;
      }
    }
    for (int i = 0; i < (int)ray_ids.size(); i++)
    {
      int ray_id = ray_ids[i];
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++) 
      {
        (*neighbour_indices)(j, ray_id) = ray_ids[indices(j, i)];
      }
    }
  }
  if (centroids || normals || dimensions || mats)
  {
    for (int i = 0; i < (int)ray_ids.size(); i++)
    {
      int ray_id = ray_ids[i];
      Eigen::Vector3d centroid;
      int num_neighbours;
      for (num_neighbours = 0; num_neighbours < search_size && indices(num_neighbours, i) != Nabo::NNSearchD::InvalidIndex; num_neighbours++){}

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(3);
      eigenSolve(ray_ids, indices, i, num_neighbours, eigen_solver, centroid);
      if (reject_back_facing_rays)
      {
        Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
        if ((ends[ray_id] - starts[ray_id]).dot(normal) > 0.0)
          normal = -normal;
        bool changed = false;
        for (int j = num_neighbours - 1; j >= 0; j--)
        {
          int id = ray_ids[indices(j, i)];
          if ((ends[id] - starts[id]).dot(normal) > 0.0)
          {
            indices(j, i) = indices(--num_neighbours, i);
            changed = true;
          }
        }
        if (changed)
        {
          eigenSolve(ray_ids, indices, i, num_neighbours, eigen_solver, centroid);
        }
      }   
      if (centroids)
        (*centroids)[ray_id] = centroid;
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
        (*dimensions)[ray_id] =
          Eigen::Vector3d(std::sqrt(eigenvals[0]), std::sqrt(eigenvals[1]), std::sqrt(eigenvals[2]));
      }
      if (mats)
        (*mats)[ray_id] = eigen_solver.eigenvectors();
    }
  }
}

// starts are required to get the normal the right way around
std::vector<Eigen::Vector3d> Cloud::generateNormals(int search_size)
{
  std::vector<Eigen::Vector3d> normals;
  getSurfels(search_size, nullptr, &normals, nullptr, nullptr, nullptr);
  return normals;
}

bool RAYLIB_EXPORT Cloud::getInfo(const std::string &file_name, Info &info)
{
  double min_s = std::numeric_limits<double>::max();
  double max_s = std::numeric_limits<double>::lowest();
  Eigen::Vector3d min_v(min_s, min_s, min_s);
  Eigen::Vector3d max_v(max_s, max_s, max_s);
  Cuboid unbounded(min_v, max_v);
  info.ends_bound = info.starts_bound = info.rays_bound = unbounded;
  info.num_rays = info.num_bounded = 0;
  info.min_time = min_s;
  info.max_time = max_s;
  info.centroid.setZero();
  info.start_pos.setZero();
  info.end_pos.setZero();
  auto find_bounds = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha > 0)
      {
        info.ends_bound.min_bound_ = minVector(info.ends_bound.min_bound_, ends[i]);
        info.ends_bound.max_bound_ = maxVector(info.ends_bound.max_bound_, ends[i]);
        info.num_bounded++;
        info.centroid += ends[i];
      }
      info.num_rays++;
      info.starts_bound.min_bound_ = minVector(info.starts_bound.min_bound_, starts[i]);
      info.starts_bound.max_bound_ = maxVector(info.starts_bound.max_bound_, starts[i]);
      info.rays_bound.min_bound_ = minVector(info.rays_bound.min_bound_, ends[i]);
      info.rays_bound.max_bound_ = maxVector(info.rays_bound.max_bound_, ends[i]);
      if (times[i] < info.min_time)
      {
        info.start_pos = starts[i];
      }
      info.min_time = std::min(info.min_time, times[i]);
      if (times[i] > info.max_time)
      {
        info.end_pos = starts[i];
      }
      info.max_time = std::max(info.max_time, times[i]);
    }
    info.rays_bound.min_bound_ = minVector(info.rays_bound.min_bound_, info.starts_bound.min_bound_);
    info.rays_bound.max_bound_ = maxVector(info.rays_bound.max_bound_, info.starts_bound.max_bound_);
  };
  bool success = readPly(file_name, true, find_bounds, 0);
  info.centroid /= static_cast<double>(info.num_bounded);
  return success;
}


double Cloud::estimatePointSpacing(const std::string &file_name, const Cuboid &bounds, int num_points)
{
  // two-iteration estimation, modelling the point distribution by the below exponent.
  // larger exponents (towards 2.5) match thick forests, lower exponents (towards 2) match smooth terrain and surfaces
  const double cloud_exponent = 2.0;  // model num_points = (cloud_width/voxel_width)^cloud_exponent

  Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  double cloud_width = pow(extent[0] * extent[1] * extent[2], 1.0 / 3.0);  // an average
  double voxel_width = cloud_width / pow((double)num_points, 1.0 / cloud_exponent);
  voxel_width *=
    5.0;  // we want to use a larger width because this process only works when the width is an overestimation
  std::cout << "initial voxel width estimate: " << voxel_width << std::endl;
  double num_voxels = 0;
  std::set<Eigen::Vector3i, Vector3iLess> test_set;

  auto estimate_size = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                           std::vector<ray::RGBA> &colours) {
    for (unsigned int i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha == 0)
        continue;

      const Eigen::Vector3d &point = ends[i];
      Eigen::Vector3i place(int(std::floor(point[0] / voxel_width)), int(std::floor(point[1] / voxel_width)),
                            int(std::floor(point[2] / voxel_width)));
      if (test_set.insert(place).second)
      {
        num_voxels++;
      }
    }
  };
  if (!readPly(file_name, true, estimate_size, 0))
    return 0;

  double points_per_voxel = (double)num_points / num_voxels;
  double width = voxel_width / pow(points_per_voxel, 1.0 / cloud_exponent);
  std::cout << "estimated point spacing: " << width << std::endl;
  return width;
}

double Cloud::estimatePointSpacing() const
{
  // two-iteration estimation, modelling the point distribution by the below exponent.
  // larger exponents (towards 2.5) match thick forests, lower exponents (towards 2) match smooth terrain and surfaces
  const double cloud_exponent = 2.0;  // model num_points = (cloud_width/voxel_width)^cloud_exponent

  Eigen::Vector3d min_bound, max_bound;
  calcBounds(&min_bound, &max_bound, kBFEnd);
  Eigen::Vector3d extent = max_bound - min_bound;
  int num_points = 0;
  for (unsigned int i = 0; i < ends.size(); i++)
    if (rayBounded(i))
      num_points++;
  double cloud_width = pow(extent[0] * extent[1] * extent[2], 1.0 / 3.0);  // an average
  double voxel_width = cloud_width / pow((double)num_points, 1.0 / cloud_exponent);
  voxel_width *=
    5.0;  // we want to use a larger width because this process only works when the width is an overestimation
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
      if (test_set.insert(place).second)
      {
        num_voxels++;
      }
    }
  }
  double points_per_voxel = (double)num_points / num_voxels;
  double width = voxel_width / pow(points_per_voxel, 1.0 / cloud_exponent);
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

void Cloud::resize(size_t size)
{
  starts.resize(size);
  ends.resize(size);
  times.resize(size);
  colours.resize(size);
}

void Cloud::reserve(size_t size)
{
  starts.reserve(size);
  ends.reserve(size);
  times.reserve(size);
  colours.reserve(size);
}

Eigen::Array<double, 22, 1> Cloud::getMoments() const
{
  Eigen::Vector3d startMean(0, 0, 0);
  Eigen::Array3d startSigma(0, 0, 0);
  Eigen::Vector3d endMean(0, 0, 0);
  Eigen::Array3d endSigma(0, 0, 0);
  double timeMean = 0.0;
  double timeSigma = 0.0;
  Eigen::Vector4d colourMean(0, 0, 0, 0);
  Eigen::Array4d colourSigma(0, 0, 0, 0);
  for (size_t i = 0; i < ends.size(); i++)
  {
    startMean += starts[i];
    endMean += ends[i];
    timeMean += times[i];
    colourMean += Eigen::Vector4d(colours[i].red, colours[i].green, colours[i].blue, colours[i].alpha) / 255.0;
  }
  startMean /= (double)ends.size();
  endMean /= (double)ends.size();
  timeMean /= (double)ends.size();
  colourMean /= (double)ends.size();
  for (size_t i = 0; i < ends.size(); i++)
  {
    Eigen::Array3d start = (starts[i] - startMean).array();
    startSigma += start * start;
    Eigen::Array3d end = (ends[i] - endMean).array();
    endSigma += end * end;
    timeSigma += ray::sqr(times[i] - timeMean);
    Eigen::Vector4d colour(colours[i].red, colours[i].green, colours[i].blue, colours[i].alpha);
    Eigen::Array4d col = (colour / 255.0 - colourMean).array();
    colourSigma += col * col;
  }
  startSigma = (startSigma / (double)ends.size()).sqrt();
  endSigma = (endSigma / (double)ends.size()).sqrt();
  timeSigma = std::sqrt(timeSigma / (double)ends.size());
  colourSigma = (colourSigma / (double)ends.size()).sqrt();

  Eigen::Array<double, 22, 1> result;
  result << startMean, startSigma, endMean, endSigma, timeMean, timeSigma, colourMean, colourSigma;
  std::cout << "stats: ";
  for (int i = 0; i < 22; i++) std::cout << ", " << result[i];
  std::cout << std::endl;
  return result;  // Note: this is used once per cloud, returning by value is not a performance issue
}

bool Cloud::read(const std::string &file_name,
                 std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                    std::vector<double> &times, std::vector<RGBA> &colours)>
                   apply)
{
  return readPly(file_name, true, apply, 0);
}

}  // namespace ray
