// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas, Tom Lowe
#include "raytransientfilter.h"

#include "raygrid.h"
#include "rayindexingoctreenode.h"
#include "rayprogress.h"
#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB

#include <algorithm>
#include <iomanip>
#include <iostream>

// #pragma GCC optimize("O0")

using namespace ray;

void EllipsoidMark::reset(size_t id)
{
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  id_ = id;
  first_intersection_time_ = std::numeric_limits<double>::max();
  last_intersection_time_ = std::numeric_limits<double>::lowest();
  hits_ = 0u;
}

void EllipsoidMark::sortPassThroughIds()
{
  std::sort(pass_through_ids_.begin(), pass_through_ids_.end());
}

void EllipsoidMark::hit(size_t ray_id, double time)
{
  RAYLIB_UNUSED(ray_id);
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  ++hits_;
  first_intersection_time_ = std::min(time, first_intersection_time_);
  last_intersection_time_ = std::max(time, last_intersection_time_);
}


void EllipsoidMark::passThrough(size_t ray_id)
{
#if RAYLIB_WITH_TBB
  Mutex::scoped_lock guard(*lock_);
#endif  // RAYLIB_WITH_TBB
  pass_through_ids_.emplace_back(ray_id);
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
  Eigen::Matrix4d inverse_sphere_transform;
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    EllipsoidMark ellipsoid_mark(i);
#if RAYLIB_ELLIPSOID_TRANSFORM
    const Ellipsoid &ellipsoid = ellipsoids_[i];
    // Create a transfomation matrix into an ellipsoid space where it becomes a unit sphere.
    inverse_sphere_transform = Eigen::Matrix4d::Identity();
    // We use an pessimistic padding of equal to half the diagonal length of a voxel.
    const Eigen::Vector3d voxel_diag =
      0.5 * Eigen::Vector3d(config_.voxel_size, config_.voxel_size, config_.voxel_size);
    const double padding = voxel_diag.norm();
    for (int i = 0; i < 3; ++i)
    {
      Eigen::Vector3d v = ellipsoid.vectors[i];
      // Padd the vector by the voxel size. This allows use to test voxel centres against a sphere.
      double len = v.norm();
      if (len > 1e-8)
      {
        v = (v / len) * (len + padding);
      }
      inverse_sphere_transform.col(i) = v.homogeneous();
    }
    inverse_sphere_transform.col(3) = ellipsoid.pos.homogeneous();
    inverse_sphere_transform = inverse_sphere_transform.inverse();
    ellipsoid_mark.setInverseTransform(inverse_sphere_transform);
#endif  // RAYLIB_ELLIPSOID_TRANSFORM
    ellipsoids_marks_.emplace_back(ellipsoid_mark);
    progress->increment();
  }

  std::unique_ptr<IndexingOctreeNode> octree(generateEllipsoidOctree(*progress));

  // // Build a grid of ellipsoid overlaps.
  // Grid<size_t> ellipsoid_grid(box_min, box_max, config_.voxel_size);
  // generateEllipsoidGrid(ellipsoid_grid, *progress);

  // transient_marks_.resize(cloud.rayCount(), false);
  // markIntersectedEllipsoids(cloud, ellipsoid_grid, true, *progress);

  transient_marks_.resize(cloud.rayCount(), false);
  markIntersectedEllipsoids(cloud, *octree, true, *progress);

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

#if RAYLIB_ELLIPSOID_TRANSFORM
  size_t total_voxels = 0;
  size_t saved_voxels = 0;
#endif  // RAYLIB_ELLIPSOID_TRANSFORM
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    const Ellipsoid &ellipsoid = ellipsoids_[i];
#if RAYLIB_ELLIPSOID_TRANSFORM
    const EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[i];
#endif  // RAYLIB_ELLIPSOID_TRANSFORM

    Eigen::Vector3d ellipsoid_min = ellipsoid.pos - ellipsoid.extents;
    Eigen::Vector3d ellipsoid_max = ellipsoid.pos + ellipsoid.extents;

    // Add the ellipsoid to the appropriate grid cells.
    Eigen::Vector3i index_min = grid.index(ellipsoid_min, true);
    Eigen::Vector3i index_max = grid.index(ellipsoid_max, true);

    // Add to the overlapping voxels.
    // Note: this is an overestimated overlap. With RAYLIB_ELLIPSOID_TRANSFORM we try reduce the number of voxels
    // added. The ellipsoid_mark has a transformation matrix which takes the voxel centres into a space where the
    // ellipsoid is spherical. We can then check for intersection between the sphere/voxel with a simple range test.
    for (int z = index_min.z(); z <= index_max.z(); ++z)
    {
      for (int y = index_min.y(); y <= index_max.y(); ++y)
      {
        for (int x = index_min.x(); x <= index_max.x(); ++x)
        {
#if RAYLIB_ELLIPSOID_TRANSFORM
          ++total_voxels;
          // Transform the voxel into the elliposid's space where the ellipsoid is a (padded) sphere. We only
          // add to the grid if the voxel centre falls within the sphere. We've already calculated the transformation
          // matrix in ellipsoid_mark.
          Eigen::Vector3d voxel_centre_transformed =
            (ellipsoid_mark.inverseTransform() * grid.voxelCentre(x, y, z).homogeneous()).head<3>();
          if (voxel_centre_transformed.squaredNorm() <= 1.0)
#endif  // RAYLIB_ELLIPSOID_TRANSFORM
          {
            grid.insert(x, y, z, i);
          }
#if RAYLIB_ELLIPSOID_TRANSFORM
          else
          {
            ++saved_voxels;
          }
#endif  // RAYLIB_ELLIPSOID_TRANSFORM
        }
      }
    }

    progress.increment();
  }

#if RAYLIB_ELLIPSOID_TRANSFORM
  std::cout << "\nsaved / total : " << saved_voxels << '/' << total_voxels << " : "
            << 100.0 * double(saved_voxels) / double(total_voxels) << '%' << std::endl;
#endif  // RAYLIB_ELLIPSOID_TRANSFORM
}

IndexingOctreeNode *TransientFilter::generateEllipsoidOctree(Progress &progress)
{
  progress.reset("ellipsoids-octree", ellipsoids_.size());

  if (ellipsoids_.empty())
  {
    // Nothing to build.
    return nullptr;
  }

  const double voxel_size = config_.voxel_size;

  // Start with at least a 1m voxel around the first ellipsoid.
  Eigen::Vector3d bounds_min = ellipsoids_[0].pos - 0.5 * ellipsoids_[0].extents;
  Eigen::Vector3d bounds_max = ellipsoids_[0].pos + 0.5 * ellipsoids_[0].extents;

  bounds_min = bounds_min.array().floor();
  bounds_max = bounds_max.array().ceil();

  IndexingOctreeNode *root = new IndexingOctreeNode(bounds_min, bounds_max, true);

  // Build the trueaddAabb
  IndexingOctreeNode *target = nullptr;
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    bounds_min = ellipsoids_[i].pos - 0.5 * ellipsoids_[i].extents;
    bounds_max = ellipsoids_[i].pos + 0.5 * ellipsoids_[i].extents;
    // Check for growth first.
    if (!root->contains(bounds_min) || !root->contains(bounds_max))
    {
      root = static_cast<IndexingOctreeNode *>(root->expand(bounds_min, bounds_max));
    }

    target = static_cast<IndexingOctreeNode *>(root->addAabb(bounds_min, bounds_max, voxel_size));
    // TODO: assert(target)
    target->addDatum(i);
    progress.increment();
  }

  return root;
}


void TransientFilter::markIntersectedEllipsoids(const Cloud &cloud, Grid<size_t> &grid, bool self_transient,
                                                Progress &progress)
{
  progress.reset("transient-mark-ellipsoids", cloud.ends.size());

  // We walk each ray along the voxel `grid`.

  // Trace rays through the ellipsoids.
#if RAYLIB_WITH_TBB
  auto tbb_process_ray = [this, &cloud, &grid, &progress](size_t ray_id)  //
  {
    walkRay(cloud, grid, ray_id);
    progress.increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ray);
#else   // RAYLIB_WITH_TBB
  for (size_t ray_id = 0; ray_id < cloud.rayCount(); ++ray_id)
  {
    walkRay(cloud, grid, ray_id);
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


void TransientFilter::markIntersectedEllipsoids(const Cloud &cloud, const IndexingOctreeNode &octree,
                                                bool self_transient, Progress &progress)
{
  progress.reset("transient-mark-ellipsoids", cloud.ends.size());

  // We trace each ray into the octree
#if RAYLIB_WITH_TBB
  auto tbb_process_ray = [this, &cloud, &octree, &progress](size_t ray_id)  //
  {
    OctreeNode::Ray ray(cloud.starts[ray_id], cloud.ends[ray_id]);
    const auto ray_visit =
      [this, &ray_id, &cloud](const OctreeNode::Ray &ray, const OctreeNode *n, const OctreeNode * /*parent*/)  //
    {
      const IndexingOctreeNode &node = *static_cast<const IndexingOctreeNode *>(n);
      if (!node.indices().empty())
      {
        rayTouch(ray_id, cloud, ray, node.indices());
      }
    };

    octree.rayTrace(ray, ray_visit);
    progress.increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ray);
#else   // RAYLIB_WITH_TBB
  for (size_t ray_id = 0; ray_id < cloud.rayCount(); ++ray_id)
  {
    OctreeNode::Ray ray(cloud.starts[ray_id], cloud.ends[ray_id]);
    const auto ray_visit =
      [this, &ray_id, &cloud](const OctreeNode::Ray &ray, const OctreeNode *n, const OctreeNode * /*parent*/)  //
    {
      const IndexingOctreeNode &node = *static_cast<const IndexingOctreeNode *>(n);
      if (!node.indices().empty())
      {
        rayTouch(ray_id, cloud, ray, node.indices());
      }
    };

    octree.rayTrace(ray, ray_visit);
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
        ellipsoid_mark.passThrough(ray_id);
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


void TransientFilter::rayTouch(size_t ray_id, const Cloud &cloud, const OctreeNode::Ray &ray,
                               const std::vector<size_t> &ellipsoid_ids)
{
  const double pass_distance = 0.05;
  const double ratio = pass_distance * ray.length_inverse;
  const Eigen::Vector3d scaled_dir = ray.direction * ray.length;
  for (size_t ellipsoid_index : ellipsoid_ids)
  {
    const Ellipsoid &ellipsoid = ellipsoids_[ellipsoid_index];
    EllipsoidMark &ellipsoid_mark = ellipsoids_marks_[ellipsoid_index];
    // ray-ellipsoid intersection
    const Eigen::Vector3d to_sphere = ellipsoid.pos - ray.origin;
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
      ellipsoid_mark.passThrough(ray_id);
    }
    else
    {
      ellipsoid_mark.hit(ray_id, cloud.times[ray_id]);
    }
  }
}
