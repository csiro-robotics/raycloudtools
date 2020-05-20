// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas, Tom Lowe
#include "raytransientfilter.h"

#include "raygrid.h"
#include "rayprogress.h"
#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>

using namespace ray;

TransientFilter::TransientFilter(const TransientFilterConfig &config)
  : config_(config)
{}

TransientFilter::~TransientFilter() = default;

bool TransientFilter::filter(const Cloud &cloud, Progress *progress)
{
  switch (config_.strategy)
  {
  case TransientFilterStrategy::EllipseGrid:
    return filterWithEllipseGrid(cloud, progress);
  case TransientFilterStrategy::RayGrid:
    return filterWithRayGrid(cloud, progress);
  }

  return false;
}

bool TransientFilter::filterWithEllipseGrid(const Cloud &cloud, Progress *progress)
{
  // Ensure we have a value progress pointer to update. This simplifies code below.
  Progress tracker;
  if (!progress)
  {
    progress = &tracker;
  }

  clear();

  Eigen::Vector3d box_min, box_max;
  cloud.generateEllipsoids(ellipsoids_, &box_min, &box_max, progress);

  progress->reset("initialise-marks", ellipsoids_.size());
  ellipsoids_marks_.reserve(ellipsoids_.size());
  Eigen::Matrix4d inverse_sphere_transform;
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    EllipsoidMark ellipsoid_mark(i);
    ellipsoids_marks_.emplace_back(ellipsoid_mark);
    progress->increment();
  }

  Grid<size_t> ellipse_grid(box_min, box_max, config_.voxel_size);
  generateEllipseGrid(ellipse_grid, *progress);

  markIntersectedEllipsoidsWithEllipseGrid(cloud, ellipse_grid, true, *progress);

  finaliseFilter(cloud, transient_marks_);

  return true;
}


bool TransientFilter::filterWithRayGrid(const Cloud &cloud, Progress *progress)
{
  // Ensure we have a value progress pointer to update. This simplifies code below.
  Progress tracker;
  if (!progress)
  {
    progress = &tracker;
  }

  clear();

  Eigen::Vector3d bounds_min, bounds_max;
  if (!cloud.calcBounds(&bounds_min, &bounds_max, kBFEnd | kBFStart, progress))
  {
    // To calculate bounds.
    return false;
  }

  progress->reset("build-ray-grid", cloud.rayCount());
  Grid<size_t> ray_grid(bounds_min, bounds_max, config_.voxel_size);
  buildRayGrid(ray_grid, cloud.starts, cloud.ends, progress);

  Eigen::Vector3d box_min, box_max;
  cloud.generateEllipsoids(ellipsoids_, &box_min, &box_max, progress);

  // Atomic do not support assignment and construction so we can't really retain the vector memory.
  std::vector<std::atomic_bool> transient_marks(cloud.rayCount());
  markIntersectedEllipsoidsWithRayGrid(cloud, ray_grid, transient_marks, true, *progress);

  finaliseFilter(cloud, transient_marks);

  return true;
}


void TransientFilter::clear()
{
  transient_.clear();
  fixed_.clear();
  ellipsoids_.clear();
  ellipsoids_marks_.clear();
  transient_marks_.clear();
}


void TransientFilter::generateEllipseGrid(Grid<size_t> &ellipse_grid, Progress &progress)
{
  progress.reset("transient-ellipsoid-grid", ellipsoids_.size());

  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    const Ellipsoid &ellipsoid = ellipsoids_[i];

    Eigen::Vector3d ellipsoid_min = ellipsoid.pos - ellipsoid.extents;
    Eigen::Vector3d ellipsoid_max = ellipsoid.pos + ellipsoid.extents;

    // Add the ellipsoid to the appropriate grid cells.
    Eigen::Vector3i index_min = ellipse_grid.index(ellipsoid_min, true);
    Eigen::Vector3i index_max = ellipse_grid.index(ellipsoid_max, true);

    // Add to the overlapping voxels.
    // Note: this is an overestimated overlap. We hvae tried reduce the number of voxels by using an inverse
    // transfomration of the ellipsoid. This takes the voxel centres into a space where the ellipsoid is spherical. We
    // can then check for intersection between the sphere/voxel with a simple range test.
    // However, this did not yield any performance benefits and used more memory to store the matrix.
    for (int z = index_min.z(); z <= index_max.z(); ++z)
    {
      for (int y = index_min.y(); y <= index_max.y(); ++y)
      {
        for (int x = index_min.x(); x <= index_max.x(); ++x)
        {
          ellipse_grid.insert(x, y, z, i);
        }
      }
    }

    progress.increment();
  }
}


void TransientFilter::markIntersectedEllipsoidsWithEllipseGrid(const Cloud &cloud, Grid<size_t> &ellipse_grid,
                                                               bool self_transient, Progress &progress)
{
  progress.reset("transient-mark-ellipsoids", cloud.ends.size());

  transient_marks_.resize(cloud.rayCount(), false);


  // We walk each ray along the voxel `ellipse_grid`.

  // Trace rays through the ellipsoids.
#if RAYLIB_WITH_TBB
  auto tbb_process_ray = [this, &cloud, &ellipse_grid, &progress](size_t ray_id)  //
  {
    walkRay(cloud, ellipse_grid, ray_id);
    progress.increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ray);
#else   // RAYLIB_WITH_TBB
  for (size_t ray_id = 0; ray_id < cloud.rayCount(); ++ray_id)
  {
    walkRay(cloud, ellipse_grid, ray_id);
    progress.increment();
  }
#endif  // RAYLIB_WITH_TBB

  progress.reset("transient-update-ellipsoids", ellipsoids_.size());

  // Process ray results.
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    size_t num_before = 0, num_after = 0;
    Ellipsoid &ellipsoid = ellipsoids_[i];
    EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[i];

    // Increment progress immediately due to the number of continue statements below.
    progress.increment();

    // May have run multi-threaded. Sort the pass_through_ids.
    ellipsoid_mark.sortPassThroughIds();

    size_t num_rays = ellipsoid_mark.hits() + ellipsoid_mark.passThroughIds().size();
    if (num_rays == 0 || self_transient)
    {
      ellipsoid.opacity = (double)ellipsoid_mark.hits() /
                          ((double)ellipsoid_mark.hits() + (double)ellipsoid_mark.passThroughIds().size());
    }

    if (num_rays == 0 || ellipsoid.opacity == 0 || config_.num_rays_filter_threshold == 0)
    {
      continue;
    }

    if (self_transient)
    {
      ellipsoid.num_gone = ellipsoid_mark.passThroughIds().size();
      // now get some density stats...
      double misses = 0;
      for (auto &ray_id : ellipsoid_mark.passThroughIds())
      {
        if (cloud.times[ray_id] > ellipsoid_mark.lastIntersectionTime())
        {
          num_after++;
        }
        else if (cloud.times[ray_id] < ellipsoid_mark.firstIntersectionTime())
        {
          num_before++;
        }
        else
        {
          misses++;
        }
      }
      double h = double(ellipsoid_mark.hits()) + 1e-8 - 1.0;  // subtracting 1 gives an unbiased opacity estimate
      ellipsoid.opacity = h / (h + misses);
      ellipsoid.num_gone = num_before + num_after;
    }
    else  // compare to other cloud
    {
      if (ellipsoid_mark.passThroughIds().size() > 0)
      {
        if (cloud.times[ellipsoid_mark.passThroughIds()[0]] > ellipsoid.time)
        {
          num_after = ellipsoid_mark.passThroughIds().size();
        }
        else
        {
          num_before = ellipsoid_mark.passThroughIds().size();
        }
      }
    }

    double sequence_length = config_.num_rays_filter_threshold / ellipsoid.opacity;
    bool remove_ellipsoid = false;
    if (config_.merge_type == MergeType::Oldest || config_.merge_type == MergeType::Newest)
    {
      if (double(std::max(num_before, num_after)) < sequence_length)
      {
        continue;
      }

      if (config_.merge_type == MergeType::Oldest)
      {
        // if false then remove numAfter rays if > seqLength
        remove_ellipsoid = double(num_before) >= sequence_length;
      }
      else if (config_.merge_type == MergeType::Newest)
      {
        // if false then remove numBefore rays if > seqLength
        remove_ellipsoid = double(num_after) >= sequence_length;
      }
    }
    else
    {
      // we use sum rather than max below, because it better picks out moving objects that may have some
      // pass through rays before and after the hit points.
      if (double(num_before + num_after) < sequence_length)
      {
        continue;
      }
      // min is remove ellipsoid, max is remove ray
      remove_ellipsoid = config_.merge_type == MergeType::Mininum;
    }

    if (remove_ellipsoid)
    {
      ellipsoid.transient = true;
    }
    // if we don't remove the ellipsoid then we should remove numBefore and numAfter rays if they're greater than
    // sequence length
    else
    {
      double d = 0.0;
      for (size_t j = 0; j < ellipsoid_mark.passThroughIds().size(); j++)
      {
        d += ellipsoid.opacity;
        if (d >= 1.0)
        {
          d--;
        }
        else
        {
          continue;
        }

        size_t i = ellipsoid_mark.passThroughIds()[j];
        if (!self_transient || cloud.times[i] < ellipsoid_mark.firstIntersectionTime() ||
            cloud.times[i] > ellipsoid_mark.lastIntersectionTime())
        {
          // remove ray i
          transient_marks_[i] = true;
        }
      }
    }
  };
}


void TransientFilter::markIntersectedEllipsoidsWithRayGrid(const Cloud &cloud, Grid<size_t> &ray_grid,
                                                           std::vector<std::atomic_bool> &transient_marks,
                                                           bool self_transient, Progress &progress)
{
  progress.reset("transient-mark-ellipsoids", cloud.ends.size());

  // transient_marks_.resize(cloud.rayCount());
  // for (size_t i = 0; i < cloud.rayCount(); ++i)
  // {
  //   transient_marks_.emplace_back();
  // }

  // Check each ellipsoid against the ray grid for intersections.
#if RAYLIB_WITH_TBB
  auto tbb_process_ellipsoid =
    [this, &cloud, &ray_grid, &transient_marks, &progress, self_transient](size_t ellipsoid_id)  //
  {
    checkEllipsoid(&ellipsoids_[ellipsoid_id], &transient_marks, cloud, ray_grid, config_.num_rays_filter_threshold,
                   config_.merge_type, self_transient);
    progress.increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ellipsoid);
#else   // RAYLIB_WITH_TBB
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    checkEllipsoid(&ellipsoids_[i], &transient_marks, cloud, ray_grid, config_.num_rays_filter_threshold,
                   config_.merge_type, self_transient);
    progress.increment();
  }
#endif  // RAYLIB_WITH_TBB
}


template <typename BOOL>
void TransientFilter::finaliseFilter(const Cloud &cloud, const std::vector<BOOL> &transient_marks)
{
  // Lastly, generate the new ray clouds from this sphere information
  for (size_t i = 0; i < ellipsoids_.size(); i++)
  {
    RGBA col = cloud.colours[i];
    if (config_.colour_cloud)
    {
      col.red = (uint8_t)((1.0 - ellipsoids_[i].planarity) * 255.0);
      col.blue = (uint8_t)(ellipsoids_[i].opacity * 255.0);
      col.green = (uint8_t)((double)ellipsoids_[i].num_gone / ((double)ellipsoids_[i].num_gone + 10.0) * 255.0);
    }

    if (ellipsoids_[i].transient || transient_marks[i])
    {
      transient_.starts.emplace_back(cloud.starts[i]);
      transient_.ends.emplace_back(cloud.ends[i]);
      transient_.times.emplace_back(cloud.times[i]);
      transient_.colours.emplace_back(col);
    }
    else
    {
      fixed_.starts.emplace_back(cloud.starts[i]);
      fixed_.ends.emplace_back(cloud.ends[i]);
      fixed_.times.emplace_back(cloud.times[i]);
      fixed_.colours.emplace_back(col);
    }
  }
}


void TransientFilter::walkRay(const Cloud &cloud, const Grid<size_t> &ellipse_grid, size_t ray_id)
{
  const auto visit = [this, &ellipse_grid, &cloud, ray_id]  //
    (const Grid<size_t> &, const Eigen::Vector3i &cell_index,
     const GridRayInfo &ray_info)  //
  {
    RAYLIB_UNUSED(ray_info);
    auto &cell = ellipse_grid.cell(cell_index);
    for (size_t ellipsoid_index : cell.data)
    {
      const Ellipsoid &ellipsoid = ellipsoids_[ellipsoid_index];
      EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[ellipsoid_index];

      switch (ellipsoid.intersect(cloud.starts[ray_id], cloud.ends[ray_id]))
      {
      default:
      case IntersectResult::Miss:
        // misses the ellipsoid
        break;
      case IntersectResult::Passthrough:
        ellipsoid_mark.passThrough(ray_id);
        break;
      case IntersectResult::Hit:
        ellipsoid_mark.hit(ray_id, cloud.times[ray_id]);
        break;
      }
    }
  };

  if (cloud.rayBounded(ray_id))
  {
    ellipse_grid.walkVoxels(cloud.starts[ray_id], cloud.ends[ray_id], visit, true);
  }
}


void TransientFilter::checkEllipsoid(Ellipsoid *ellipsoid, std::vector<std::atomic_bool> *transients,
                                     const Cloud &cloud, const Grid<size_t> &ray_grid, double num_rays,
                                     MergeType merge_type, bool self_transient) const
{
  if (ellipsoid->transient)
  {
    // Already marked for removal. Nothing to do.
    return;
  }

  // unbounded rays cannot be a transient object
  if (ellipsoid->extents == Eigen::Vector3d::Zero())
  {
    return;
  }

  // get all the rays that overlap this ellipsoid
  const Eigen::Vector3d ellipsoid_bounds_min =
    (ellipsoid->pos - ellipsoid->extents - ray_grid.box_min) / ray_grid.voxel_width;
  const Eigen::Vector3d ellipsoid_bounds_max =
    (ellipsoid->pos + ellipsoid->extents - ray_grid.box_min) / ray_grid.voxel_width;

  if (ellipsoid_bounds_max[0] < 0.0 || ellipsoid_bounds_max[1] < 0.0 || ellipsoid_bounds_max[2] < 0.0)
  {
    return;
  }

  if (ellipsoid_bounds_min[0] > (double)ray_grid.dims[0] - 1.0 ||
      ellipsoid_bounds_min[1] > (double)ray_grid.dims[1] - 1.0 ||
      ellipsoid_bounds_min[2] > (double)ray_grid.dims[2] - 1.0)
  {
    // Out of bounds against the ray grid.
    return;
  }

  Eigen::Vector3i bmin = maxVector(Eigen::Vector3i(0, 0, 0), Eigen::Vector3i(ellipsoid_bounds_min.cast<int>()));
  Eigen::Vector3i bmax = minVector(Eigen::Vector3i(ellipsoid_bounds_max.cast<int>()),
                                   Eigen::Vector3i(ray_grid.dims[0] - 1, ray_grid.dims[1] - 1, ray_grid.dims[2] - 1));

  std::set<size_t> ray_ids;
  for (int x = bmin[0]; x <= bmax[0]; x++)
  {
    for (int y = bmin[1]; y <= bmax[1]; y++)
    {
      for (int z = bmin[2]; z <= bmax[2]; z++)
      {
        auto &ray_list = ray_grid.cell(x, y, z).data;
        for (auto &ray_id : ray_list)
        {
          ray_ids.emplace(ray_id);
        }
      }
    }
  }

  double first_intersection_time = std::numeric_limits<double>::max();
  double last_intersection_time = std::numeric_limits<double>::lowest();
  unsigned hits = 0;
  std::vector<size_t> pass_through_ids;
  for (auto &ray_id : ray_ids)
  {
    switch (ellipsoid->intersect(cloud.starts[ray_id], cloud.ends[ray_id]))
    {
    default:
    case IntersectResult::Miss:
      // misses the ellipsoid
      break;
    case IntersectResult::Passthrough:
      pass_through_ids.push_back(ray_id);
      break;
    case IntersectResult::Hit:
      ++hits;
      first_intersection_time = std::min(first_intersection_time, cloud.times[ray_id]);
      last_intersection_time = std::max(last_intersection_time, cloud.times[ray_id]);
      break;
    }
  }

  size_t num_before = 0, num_after = 0;
  ellipsoid->num_rays = hits + pass_through_ids.size();
  if (num_rays == 0 || self_transient)
  {
    ellipsoid->opacity = (double)hits / ((double)hits + (double)pass_through_ids.size());
  }
  if (ellipsoid->num_rays == 0 || ellipsoid->opacity == 0 || num_rays == 0)
  {
    return;
  }
  if (self_transient)
  {
    ellipsoid->num_gone = pass_through_ids.size();
    // now get some density stats...
    double misses = 0;
    for (auto &ray_id : pass_through_ids)
    {
      if (cloud.times[ray_id] > last_intersection_time)
      {
        num_after++;
      }
      else if (cloud.times[ray_id] < first_intersection_time)
      {
        num_before++;
      }
      else
      {
        misses++;
      }
    }
    double h = hits + 1e-8 - 1.0;  // subtracting 1 gives an unbiased opacity estimate
    ellipsoid->opacity = h / (h + misses);
    ellipsoid->num_gone = num_before + num_after;
  }
  else  // compare to other cloud
  {
    if (pass_through_ids.size() > 0)
    {
      if (cloud.times[pass_through_ids[0]] > ellipsoid->time)
      {
        num_after = pass_through_ids.size();
      }
      else
      {
        num_before = pass_through_ids.size();
      }
    }
  }

  double sequence_length = num_rays / ellipsoid->opacity;
  int remove_ellipsoid = false;
  // type:
  // - 0 oldest
  // - 1 newest
  // - 2 min
  // - 3 max
  if (merge_type == MergeType::Oldest || merge_type == MergeType::Newest)
  {
    if (double(std::max(num_before, num_after)) < sequence_length)
    {
      return;
    }
    if (merge_type == MergeType::Oldest)
    {  // if false then remove numAfter rays if > seqLength
      remove_ellipsoid = double(num_before) >= sequence_length;
    }
    else  // Newest
    {
      // if false then remove numBefore rays if > seqLength
      remove_ellipsoid = double(num_after) >= sequence_length;
    }
  }
  else  // min/max
  {
    // we use sum rather than max below, because it better picks out moving objects that may have some
    // pass through rays before and after the hit points.
    if (double(num_before + num_after) < sequence_length)
    {
      return;
    }
    // min is remove ellipsoid, max is remove ray
    remove_ellipsoid = merge_type == MergeType::Maximum;
  }

  if (remove_ellipsoid)
  {
    ellipsoid->transient = true;
  }
  else
  {
    // if we don't remove the ellipsoid then we should remove numBefore and numAfter rays if they're greater than
    // sequence length
    double d = 0.0;
    for (size_t j = 0; j < pass_through_ids.size(); ++j)
    {
      d += ellipsoid->opacity;
      if (d >= 1.0)
      {
        d--;
      }
      else
      {
        continue;
      }

      size_t ray_id = pass_through_ids[j];
      if (!self_transient || cloud.times[ray_id] < first_intersection_time ||
          cloud.times[ray_id] > last_intersection_time)
      {
        // remove ray i
        (*transients)[ray_id] = true;
      }
    }
  }
}
