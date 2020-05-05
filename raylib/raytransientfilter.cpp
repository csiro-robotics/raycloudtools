// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas, Tom Lowe
#include "raytransientfilter.h"

#if TRANSIENT_FILTER

#include "raygrid.h"
#include "rayprogress.h"
#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <tbb/parallel_for.h>
#endif // RAYLIB_WITH_TBB

#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace ray;

void EllipsoidMark::hit(size_t ray_id, double time)
{
  RAYLIB_UNUSED(ray_id);
#if RAYLIB_WITH_TBB
  std::unique_lock<tbb::mutex> guard(*lock);
#endif // RAYLIB_WITH_TBB
  ++hits;
  first_intersection_time = std::min(time, first_intersection_time);
  last_intersection_time = std::max(time, last_intersection_time);
}


void EllipsoidMark::passthrough(size_t ray_id)
{
#if RAYLIB_WITH_TBB
  std::unique_lock<tbb::mutex> guard(*lock);
#endif // RAYLIB_WITH_TBB
  pass_through_ids.emplace_back(ray_id);
}


TransientFilter::TransientFilter(const TransientFilterConfig &config)
  : config_(config)
{}

TransientFilter::~TransientFilter() = default;

void TransientFilter::filter(const Cloud &cloud, Progress *progress)
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

  progress->reset("build-markers", ellipsoids_.size());

  ellipsoids_marks_.reserve(ellipsoids_.size());
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    EllipsoidMark ellipsoid_mark(i);
    ellipsoids_marks_.emplace_back(ellipsoid_mark);
    progress->increment();
  }

  // Build a grid of ellipsoid overlaps.
  Grid<size_t> ellipsoid_grid(box_min, box_max, config_.voxel_size);
  generateEllipsoidGrid(ellipsoid_grid, *progress);

  transient_marks_.resize(cloud.rayCount(), false);
  markIntersectedEllipsoids(cloud, ellipsoid_grid, true, *progress);

  finaliseFilter(cloud);
}


void TransientFilter::clear()
{
  transient_.clear();
  fixed_.clear();
  ellipsoids_.clear();
  ellipsoids_marks_.clear();
  transient_marks_.clear();
}


void TransientFilter::generateEllipsoidGrid(Grid<size_t> &grid, Progress &progress)
{
  progress.reset("transient-ellipsoid-grid", ellipsoids_.size());

  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    const Ellipsoid &ellipsoid = ellipsoids_[i];

    Eigen::Vector3d ellipsoid_min = ellipsoid.pos - ellipsoid.extents;
    Eigen::Vector3d ellipsoid_max = ellipsoid.pos + ellipsoid.extents;

    // Add the ellipsoid to the appropriate grid cells.
    Eigen::Vector3i index_min = grid.index(ellipsoid_min, true);
    Eigen::Vector3i index_max = grid.index(ellipsoid_max, true);

    // Add to the overlapping voxels.
    // Note: this is an overestimated overlap.
    // TODO: (KS) There may be some utility in a more precise intersection with the ellipsoid here. One option is to do
    // the following:
    // 1. Pad the ellipsoid size by the voxel size.
    // 2. Generate a transformation matrix from cloud/ray space into the ellipsoid space. This transformation is based
    //  on:
    //    - translation : -1 * ellipsoid centre.
    //    - rotation : the inverse ellipsoid rotation matrix (from its vectors)
    //    - scale : the inverse scale of the ellipsoid (from its vectors)
    // 3. Take each voxel centre and transform it by the matrix calculated above. This essentially puts it into a space
    //  where the ellipsoid is a unit sphere at the origin.
    // 4. Reject the voxel (centre) if its transformed position is more than 1 unit away from the origin.
    for (int z = index_min.z(); z <= index_max.z(); ++z)
    {
      for (int y = index_min.y(); y <= index_max.y(); ++y)
      {
        for (int x = index_min.x(); x <= index_max.x(); ++x)
        {
          grid.insert(x, y, z, i);
        }
      }
    }

    progress.increment();
  }
}


void TransientFilter::markIntersectedEllipsoids(const Cloud &cloud, Grid<size_t> &grid, bool self_transient,
                                                Progress &progress)
{
  progress.reset("transient-mark-ellipsoids", cloud.ends.size());

  // TODO: support threading building blocks (TBB) for multi-threading.

  // We walk each ray along the voxel `grid`.

  // Trace rays through the ellipsoids.
#if RAYLIB_WITH_TBB
  auto tbb_process_ray = [this, &cloud, &grid, &progress] (size_t ray_id) // 
  {
    walkRay(cloud, grid, ray_id);
    progress.increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ray);
#else  // RAYLIB_WITH_TBB
  for (size_t ray_id = 0; ray_id < cloud.rayCount(); ++ray_id)
  {
    walkRay(cloud, grid, ray_id);
    progress.increment();
  }
#endif // RAYLIB_WITH_TBB

  progress.reset("transient-update-ellipsoids", ellipsoids_.size());

  // Process ray results.
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    size_t num_before = 0, num_after = 0;
    Ellipsoid &ellipsoid = ellipsoids_[i];
    EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[i];

    // May have run multi-threaded. Sort the pass_through_ids.
    std::sort(ellipsoid_mark.pass_through_ids.begin(), ellipsoid_mark.pass_through_ids.end());

    size_t num_rays = ellipsoid_mark.hits + ellipsoid_mark.pass_through_ids.size();
    if (num_rays == 0 || self_transient)
    {
      ellipsoid.opacity =
        (double)ellipsoid_mark.hits / ((double)ellipsoid_mark.hits + (double)ellipsoid_mark.pass_through_ids.size());
    }

    if (num_rays == 0 || ellipsoid.opacity == 0 || config_.num_rays_filter_threshold == 0)
    {
      continue;
    }

    if (self_transient)
    {
      ellipsoid.num_gone = ellipsoid_mark.pass_through_ids.size();
      // now get some density stats...
      double misses = 0;
      for (auto &ray_id : ellipsoid_mark.pass_through_ids)
      {
        if (cloud.times[ray_id] > ellipsoid_mark.last_intersection_time)
        {
          num_after++;
        }
        else if (cloud.times[ray_id] < ellipsoid_mark.first_intersection_time)
        {
          num_before++;
        }
        else
        {
          misses++;
        }
      }
      double h = double(ellipsoid_mark.hits) + 1e-8 - 1.0;  // subtracting 1 gives an unbiased opacity estimate
      ellipsoid.opacity = h / (h + misses);
      ellipsoid.num_gone = num_before + num_after;
    }
    else  // compare to other cloud
    {
      if (ellipsoid_mark.pass_through_ids.size() > 0)
      {
        if (cloud.times[ellipsoid_mark.pass_through_ids[0]] > ellipsoid.time)
        {
          num_after = ellipsoid_mark.pass_through_ids.size();
        }
        else
        {
          num_before = ellipsoid_mark.pass_through_ids.size();
        }
      }
    }

    double sequence_length = config_.num_rays_filter_threshold / ellipsoid.opacity;
    bool remove_ellipsoid = false;
    if (config_.merge_type == TransientFilterType::Oldest || config_.merge_type == TransientFilterType::Newest)
    {
      if (double(std::max(num_before, num_after)) < sequence_length)
      {
        continue;
      }

      if (config_.merge_type == TransientFilterType::Oldest)
      {
        // if false then remove numAfter rays if > seqLength
        remove_ellipsoid = double(num_before) >= sequence_length;
      }
      else if (config_.merge_type == TransientFilterType::Newest)
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
      remove_ellipsoid = config_.merge_type == TransientFilterType::Mininum;
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
      for (size_t j = 0; j < ellipsoid_mark.pass_through_ids.size(); j++)
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

        size_t i = ellipsoid_mark.pass_through_ids[j];
        if (!self_transient || cloud.times[i] < ellipsoid_mark.first_intersection_time ||
            cloud.times[i] > ellipsoid_mark.last_intersection_time)
        {
          // remove ray i
          transient_marks_[i] = true;
        }
      }
    }
    progress.increment();
  };
}


void TransientFilter::finaliseFilter(const Cloud &cloud)
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

    if (ellipsoids_[i].transient || transient_marks_[i])
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


void TransientFilter::walkRay(const Cloud &cloud, const Grid<size_t> &grid, size_t ray_id)
{
  const auto visit = [this, &grid, &cloud, ray_id]  //
    (const Grid<size_t> &, const Eigen::Vector3i &cell_index,
     const GridRayInfo &ray_info)  //
  {
    auto &cell = grid.cell(cell_index);
    const Eigen::Vector3d scaled_dir = ray_info.ray_end - ray_info.ray_start;
    const double pass_distance = 0.05;
    const double ratio = pass_distance / ray_info.ray_length;

    for (size_t ellipsoid_index : cell.data)
    {
      const Ellipsoid &ellipsoid = ellipsoids_[ellipsoid_index];
      EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[ellipsoid_index];
      // ray-ellipsoid intersection
      const Eigen::Vector3d to_sphere = ellipsoid.pos - ray_info.ray_start;
      Eigen::Vector3d ray;
      Eigen::Vector3d to;

      for (int i = 0; i < 3; ++i)
      {
        ray[i] = scaled_dir.dot(ellipsoid.vectors[i]);
        to[i] = to_sphere.dot(ellipsoid.vectors[i]);
      }

      double ray_length = ray.norm();
      ray /= ray_length;

      double d = to.dot(ray);
      double dist2 = (to - ray * d).squaredNorm();

      if (dist2 > 1.0)  // misses the ellipsoid
      {
        continue;
      }

      double along_dist = sqrt(1.0 - dist2);
      if (ray_length < d - along_dist)  // doesn'o reach the ellipsoid
      {
        continue;
      }

      // last number requires rays to pass some way past the object
      bool pass_through = ray_length * (1.0 - ratio) > d + along_dist;
      if (pass_through)
      {
        ellipsoid_mark.passthrough(ray_id);
      }
      else
      {
        ellipsoid_mark.hit(ray_id, cloud.times[ray_id]);
      }
    }
  };

  if (cloud.rayBounded(ray_id))
  {
    grid.walkVoxels(cloud.starts[ray_id], cloud.ends[ray_id], visit, true);
  }
}

#endif  // TRANSIENT_FILTER
