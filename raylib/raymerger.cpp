// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas, Tom Lowe
#include "raymerger.h"

#include "raygrid.h"
#include "rayprogress.h"
#include "rayunused.h"

#if RAYLIB_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>

#if RAYLIB_WITH_TBB
// With threads we use std::atomic_bool for the transient marks. These are default initialised to false. No additional
// argument required
#define MARKER_BOOL_INIT
#else  // RAYLIB_WITH_TBB
// Without threads, we use bool for the transient marks. There is no default construction, so we must provide the
// initialisation argument
#define MARKER_BOOL_INIT , false
#endif  // RAYLIB_WITH_TBB

namespace ray
{
class EllipsoidTransientMarker
{
public:
  EllipsoidTransientMarker(size_t ray_count)
    : ray_tested(ray_count, false)
  {}
  EllipsoidTransientMarker(const EllipsoidTransientMarker &other)
    : ray_tested(other.ray_tested.size(), false)
  {}

  /// Test a single @p ellipsoid against the @p ray_grid and resolve whether it should be marked as traisient.
  /// The @p ellipsoid is considered transient if sufficient rays pass through or near it.
  ///
  /// @param ellipsoid The ellipsoid to check for transient marks.
  /// @param transient_ray_marks Array marking which rays from @p cloud are transient and should be removed.
  /// @param ray_grid The voxelised representation of @p cloud .
  /// @param num_rays Thresholding value indicating the number of nearby rays required to mark the ellipsoid as
  /// transient.
  /// @param merge_type The merging strategy.
  /// @param self_transient True when the @p ellipsoid was generated from @p cloud and we are looking for transient
  /// points within this cloud.
  void mark(Ellipsoid *ellipsoid, std::vector<Merger::Bool> *transient_ray_marks, const Cloud &cloud,
            const Grid<unsigned> &ray_grid, double num_rays, MergeType merge_type, bool self_transient,
            bool ellipsoid_cloud_first);

private:
  // Working memory.

  /// Tracks which rays have been tested. Sized to match incoming cloud ray count.
  std::vector<bool> ray_tested;
  /// Ids of ray to test.
  std::vector<unsigned> test_ray_ids;
  /// Ids of rays which intersect the ellipsoid with a @c IntersectResult::Passthrough result.
  std::vector<unsigned> pass_through_ids;
};

typedef Eigen::Matrix<double, 6, 1> Vector6i;
class Vector6iLess
{
public:
  inline bool operator()(const Vector6i &a, const Vector6i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    if (a[2] != b[2])
      return a[2] < b[2];
    if (a[3] != b[3])
      return a[3] < b[3];
    if (a[4] != b[4])
      return a[4] < b[4];
    return a[5] < b[5];
  }
};

// TODO: Make config value
const double test_width = 0.01;  // allows a minor variation when checking for similarity of rays

void rayLookup(const Cloud *cloud, std::set<Vector6i, Vector6iLess> &ray_lookup)
{
  for (size_t i = 0; i < cloud->rayCount(); i++)
  {
    const Eigen::Vector3d &point = cloud->ends[i];
    const Eigen::Vector3d &start = cloud->starts[i];
    Vector6i ray;
    for (int j = 0; j < 3; j++)
    {
      ray[j] = int(floor(start[j] / test_width));
      ray[3 + j] = int(floor(point[j] / test_width));
    }
    ray_lookup.insert(ray);
  }
}

void EllipsoidTransientMarker::mark(Ellipsoid *ellipsoid, std::vector<Merger::Bool> *transient_ray_marks,
                                    const Cloud &cloud, const Grid<unsigned> &ray_grid, double num_rays,
                                    MergeType merge_type, bool self_transient, bool ellipsoid_cloud_first)
{
  if (ellipsoid->transient)
  {
    // Already marked for removal. Nothing to do.
    return;
  }

  // unbounded rays cannot be a transient object
  if (ellipsoid->extents == Eigen::Vector3f::Zero())
  {
    return;
  }

  if (ray_tested.size() != cloud.rayCount())
  {
    ray_tested.resize(cloud.rayCount(), false);
  }
  test_ray_ids.clear();
  pass_through_ids.clear();

  // get all the rays that overlap this ellipsoid
  const Eigen::Vector3d ellipsoid_bounds_min =
    (ellipsoid->pos - ellipsoid->extents.cast<double>() - ray_grid.box_min) / ray_grid.voxel_width;
  const Eigen::Vector3d ellipsoid_bounds_max =
    (ellipsoid->pos + ellipsoid->extents.cast<double>() - ray_grid.box_min) / ray_grid.voxel_width;

  if (ellipsoid_bounds_max[0] < 0.0 || ellipsoid_bounds_max[1] < 0.0 || ellipsoid_bounds_max[2] < 0.0)
  {
    return;
  }

  if (ellipsoid_bounds_min[0] >= (double)ray_grid.dims[0] || ellipsoid_bounds_min[1] >= (double)ray_grid.dims[1] ||
      ellipsoid_bounds_min[2] >= (double)ray_grid.dims[2])
  {
    // Out of bounds against the ray grid.
    return;
  }

  Eigen::Vector3i bmin = maxVector(Eigen::Vector3i(0, 0, 0), Eigen::Vector3i(ellipsoid_bounds_min.cast<int>()));
  Eigen::Vector3i bmax = minVector(Eigen::Vector3i(ellipsoid_bounds_max.cast<int>()),
                                   Eigen::Vector3i(ray_grid.dims[0] - 1, ray_grid.dims[1] - 1, ray_grid.dims[2] - 1));

  for (int x = bmin[0]; x <= bmax[0]; x++)
  {
    for (int y = bmin[1]; y <= bmax[1]; y++)
    {
      for (int z = bmin[2]; z <= bmax[2]; z++)
      {
        auto &ray_list = ray_grid.cell(x, y, z).data;
        for (auto &ray_id : ray_list)
        {
          if (ray_tested[ray_id])
          {
            continue;
          }
          ray_tested[ray_id] = true;
          test_ray_ids.push_back(ray_id);
        }
      }
    }
  }

  double first_intersection_time = std::numeric_limits<double>::max();
  double last_intersection_time = std::numeric_limits<double>::lowest();
  unsigned hits = 0;
  for (auto &ray_id : test_ray_ids)
  {
    ray_tested[ray_id] = false;

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
    ellipsoid->opacity = (float)hits / ((float)hits + (float)pass_through_ids.size());
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
    ellipsoid->opacity = static_cast<float>(h / (h + misses));
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
      // TODO: even a tiny bit of translucency will make a single ray not enough
      return;
    }
    if (merge_type == MergeType::Order)  // order
    {
      // if the cloud containing the ellipsoid is the first cloud in order, then remove the ray
      remove_ellipsoid = !ellipsoid_cloud_first;
    }
    else
    {
      // min is remove ellipsoid, max is remove ray
      remove_ellipsoid = merge_type == MergeType::Mininum;
    }
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

      unsigned ray_id = pass_through_ids[j];
      if (!self_transient || cloud.times[ray_id] < first_intersection_time ||
          cloud.times[ray_id] > last_intersection_time)
      {
        // remove ray i
        (*transient_ray_marks)[ray_id] = true;
      }
    }
  }
}

Merger::Merger(const MergerConfig &config)
  : config_(config)
{}

Merger::~Merger() = default;

bool Merger::filter(const Cloud &cloud, Progress *progress)
{
  // Ensure we have a value progress pointer to update. This simplifies code below.
  Progress tracker;
  if (!progress)
  {
    progress = &tracker;
  }

  clear();

  Eigen::Vector3d bounds_min, bounds_max;
  generateEllipsoids(&ellipsoids_, &bounds_min, &bounds_max, cloud, progress);

  const double voxel_size = voxelSizeForCloud(cloud);
  if (config_.voxel_size == 0)
  {
    std::cout << "estimated required voxel size: " << voxel_size << std::endl;
  }

  Grid<unsigned> ray_grid(bounds_min, bounds_max, voxel_size);
  seedRayGrid(&ray_grid, cloud);
  fillRayGrid(&ray_grid, cloud, progress);

  // Atomic do not support assignment and construction so we can't really retain the vector memory.
  std::vector<Bool> transient_ray_marks(cloud.rayCount() MARKER_BOOL_INIT);
  markIntersectedEllipsoids(cloud, ray_grid, &transient_ray_marks, config_.num_rays_filter_threshold, true, progress);

  finaliseFilter(cloud, transient_ray_marks);

  progress->end();

  return true;
}

bool Merger::mergeMultiple(std::vector<Cloud> &clouds, Progress *progress)
{
  // Ensure we have a value progress pointer to update. This simplifies code below.
  Progress tracker;
  if (!progress)
  {
    progress = &tracker;
  }

  clear();

  std::vector<Grid<unsigned>> grids(clouds.size());
  for (size_t c = 0; c < clouds.size(); c++)
  {
    const double voxel_size = voxelSizeForCloud(clouds[c]);
    if (config_.voxel_size == 0)
    {
      std::cout << "estimated required voxel size for cloud " << c << ": " << voxel_size << std::endl;
    }
    grids[c].init(clouds[c].calcMinBound(), clouds[c].calcMaxBound(), voxel_size);
    for (size_t d = 0; d < clouds.size(); d++)
    {
      seedRayGrid(&grids[c], clouds[d]);
    }
  }

  for (size_t c = 0; c < clouds.size(); c++)
  {
    fillRayGrid(&grids[c], clouds[c], progress);
  }  

  std::vector<std::vector<Bool>> transient_ray_marks;
  transient_ray_marks.reserve(clouds.size());
  for (size_t c = 0; c < clouds.size(); c++)
  {
    transient_ray_marks.emplace_back(std::vector<Bool>(clouds[c].rayCount() MARKER_BOOL_INIT));
  }

  // now for each cloud, look for other clouds that penetrate it
  for (size_t c = 0; c < clouds.size(); c++)
  {
    generateEllipsoids(&ellipsoids_, nullptr, nullptr, clouds[c], progress);
    // just set opacity
    markIntersectedEllipsoids(clouds[c], grids[c], &transient_ray_marks[c], 0, false, progress);

    for (size_t d = 0; d < clouds.size(); d++)
    {
      if (d == c)
      {
        continue;
      }
      const bool ellipsoid_cloud_first = c < d;  // used when argument order of the files is the merge type
      // use ellipsoid opacity to set transient flag true on transients
      markIntersectedEllipsoids(clouds[d], grids[d], &transient_ray_marks[d], config_.num_rays_filter_threshold, false,
                                progress, ellipsoid_cloud_first);
    }

    for (size_t i = 0; i < clouds[c].rayCount(); i++)
    {
      if (ellipsoids_[i].transient)
      {
        transient_ray_marks[c][i] = true;
      }
    }
  }

  for (size_t c = 0; c < clouds.size(); c++)
  {
    auto &cloud = clouds[c];
    for (size_t i = 0; i < cloud.rayCount(); i++)
    {
      if (transient_ray_marks[c][i])
      {
        difference_.addRay(cloud, i);
      }
      else
      {
        fixed_.addRay(cloud, i);
      }
    }
  }

  return true;
}

bool Merger::mergeThreeWay(const Cloud &base_cloud, Cloud &cloud1, Cloud &cloud2, Progress *progress)
{
  // The 3-way merge is similar to those performed on text files for version control systems. It attempts to apply the
  // changes in both cloud 1 and cloud2 (compared to base_cloud). When there is a conflict (different changes in the
  // same location) it resolves that according to the selected merge_type. unlike with text, a change requires a small
  // threshold, since positions are floating point values. In our case, we define a ray as unchanged when the start and
  // end points are within the same small voxel as they were in base_cloud. so the threshold is test_width.

  // generate quick lookup for the existance of a particular (quantised) ray
  Cloud *clouds[2] = { &cloud1, &cloud2 };
  std::set<Vector6i, Vector6iLess> base_ray_lookup;
  rayLookup(&base_cloud, base_ray_lookup);
  std::set<Vector6i, Vector6iLess> ray_lookups[2];
  for (int c = 0; c < 2; c++) rayLookup(clouds[c], ray_lookups[c]);

  std::cout << "set size " << ray_lookups[0].size() << ", " << ray_lookups[1].size() << ", " << base_ray_lookup.size()
            << std::endl;

  // now remove all similar rays to base_cloud and put them in the final cloud:
  int preferred_cloud = clouds[0]->times[0] > clouds[1]->times[0] ? 0 : 1;
  size_t u = 0;
  for (int c = 0; c < 2; c++)
  {
    Cloud &cloud = *clouds[c];
    for (size_t i = 0; i < cloud.rayCount(); i++)
    {
      Eigen::Vector3d &point = cloud.ends[i];
      Eigen::Vector3d &start = cloud.starts[i];
      Vector6i ray;
      for (int j = 0; j < 3; j++)
      {
        ray[j] = int(std::floor(start[j] / test_width));
        ray[3 + j] = int(std::floor(point[j] / test_width));
      }
      int other = 1 - c;
      // if the ray is in cloud1 and cloud2 there is no contention, so add the ray to the result
      if (ray_lookups[other].find(ray) != ray_lookups[other].end())
      {
        if (c == preferred_cloud)
        {
          fixed_.addRay(start, point, cloud.times[i], cloud.colours[i]);
          u++;
        }
      }
      // we want to run the combine (which revolves conflicts) on only the changed parts
      // so we want to keep only the changes for cloud[0] and cloud[1]...
      // which means removing rays that aren't changed:
      if (base_ray_lookup.find(ray) != base_ray_lookup.end())
      {
        cloud.starts[i] = cloud.starts.back();
        cloud.starts.pop_back();
        cloud.ends[i] = cloud.ends.back();
        cloud.ends.pop_back();
        cloud.times[i] = cloud.times.back();
        cloud.times.pop_back();
        cloud.colours[i] = cloud.colours.back();
        cloud.colours.pop_back();
        i--;
      }
    }
  }
  std::cout << u << " unaltered rays have been moved into combined cloud" << std::endl;
  std::cout << clouds[0]->rayCount() << " and " << clouds[1]->rayCount() << " rays to combine, that are different"
            << std::endl;
#if defined VERBOSE_MERGE
  fixed_.save("common_rays.ply");
  clouds[0]->save("changes_0.ply");
  clouds[1]->save("changes_1.ply");
#endif
  // This keeps all data only where there are conflicts.
  if (config_.merge_type == MergeType::All)
  {
    for (int c = 0; c < 2; c++)
    {
      fixed_.starts.insert(fixed_.starts.end(), clouds[c]->starts.begin(), clouds[c]->starts.end());
      fixed_.ends.insert(fixed_.ends.end(), clouds[c]->ends.begin(), clouds[c]->ends.end());
      fixed_.times.insert(fixed_.times.end(), clouds[c]->times.begin(), clouds[c]->times.end());
      fixed_.colours.insert(fixed_.colours.end(), clouds[c]->colours.begin(), clouds[c]->colours.end());
    }
    return true;
  }
  // otherwise we run combine on the altered clouds
  // first, grid the rays for fast lookup
  Grid<unsigned> grids[2];
  for (int c = 0; c < 2; c++)
  {
    grids[c].init(clouds[c]->calcMinBound(), clouds[c]->calcMaxBound(), voxelSizeForCloud(*clouds[c]));
    seedRayGrid(&grids[c], *clouds[0]); // to only fill rays in voxels occupied by cloud 0 or 1
    seedRayGrid(&grids[c], *clouds[1]);
    fillRayGrid(&grids[c], *clouds[c], progress);
  }

  std::vector<Bool> transients[2] = { std::vector<Bool>(clouds[0]->rayCount() MARKER_BOOL_INIT),
                                      std::vector<Bool>(clouds[1]->rayCount() MARKER_BOOL_INIT) };
  // now for each cloud, represent the end points as ellipsoids, and ray cast the other cloud's rays against it
  for (int c = 0; c < 2; c++)
  {
    if (clouds[c]->rayCount() == 0)
    {
      continue;
    }
    generateEllipsoids(&ellipsoids_, nullptr, nullptr, *clouds[c]);

    // just set opacity
    markIntersectedEllipsoids(*clouds[c], grids[c], &transients[c], 0, false, progress);

    const int d = 1 - c;
    const bool ellipsoid_cloud_first = c < d;  // used when argument order of the files is the merge type
    // use ellipsoid opacity to set transient flag true on transients (intersected ellipsoids)
    markIntersectedEllipsoids(*clouds[d], grids[d], &transients[d], config_.num_rays_filter_threshold, false, progress,
                              ellipsoid_cloud_first);

    for (size_t i = 0; i < ellipsoids_.size(); i++)
    {
      if (ellipsoids_[i].transient)
      {
        transients[c][i] = true;
      }
    }
  }
  for (int c = 0; c < 2; c++)
  {
    auto &cloud = *clouds[c];
    size_t removed_count = 0;
    for (size_t i = 0; i < transients[c].size(); i++)
    {
      if (!transients[c][i])
      {
        fixed_.addRay(cloud, i);
      }
      else
      {
        removed_count++;  // we aren't storing the differences. No current demand for this.
      }
    }
    std::cout << removed_count << " removed rays, " << fixed_.rayCount() << " fixed rays." << std::endl;
  }

  return true;
}

void Merger::clear()
{
  difference_.clear();
  fixed_.clear();
  ellipsoids_.clear();
}

void Merger::seedRayGrid(Grid<unsigned> *grid, const Cloud &cloud)
{
  const auto seed_voxels = [grid, &cloud](unsigned i)
  {
    Eigen::Vector3d end = (cloud.ends[i] - grid->box_min) / grid->voxel_width;
    Eigen::Vector3i index((int)floor(end[0]), (int)floor(end[1]), (int)floor(end[2]));
    grid->addCell(index);
  };
#if RAYLIB_PARALLEL_GRID
  tbb::parallel_for<unsigned>(0u, unsigned(cloud.rayCount()), seed_voxels);
#else   // RAYLIB_PARALLEL_GRID
  const unsigned int count = static_cast<unsigned int>(cloud.rayCount());
//  #pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < count; ++i)
  {
    seed_voxels(i);
  }
#endif  // RAYLIB_PARALLEL_GRID
}

void Merger::fillRayGrid(Grid<unsigned> *grid, const Cloud &cloud, Progress *progress)
{
  if (progress)
  {
    progress->begin("fillRayGrid", cloud.rayCount());
  }

  const auto add_ray = [grid, &cloud, progress](unsigned i)  //
  {
    Eigen::Vector3d dir = cloud.ends[i] - cloud.starts[i];
    Eigen::Vector3d dir_sign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Eigen::Vector3d start = (cloud.starts[i] - grid->box_min) / grid->voxel_width;
    Eigen::Vector3d end = (cloud.ends[i] - grid->box_min) / grid->voxel_width;
    Eigen::Vector3i start_index((int)floor(start[0]), (int)floor(start[1]), (int)floor(start[2]));
    Eigen::Vector3i end_index((int)floor(end[0]), (int)floor(end[1]), (int)floor(end[2]));
    double length_sqr = (end_index - start_index).squaredNorm();
    Eigen::Vector3i index = start_index;
    for (;;)
    {
      grid->insertIfCellExists(index, i);
      if (index == end_index || (index - start_index).squaredNorm() > length_sqr)
      {
        break;
      }
      Eigen::Vector3d mid =
        grid->box_min + grid->voxel_width * Eigen::Vector3d(index[0] + 0.5, index[1] + 0.5, index[2] + 0.5);
      Eigen::Vector3d next_boundary = mid + 0.5 * grid->voxel_width * dir_sign;
      Eigen::Vector3d delta = next_boundary - cloud.starts[i];
      Eigen::Vector3d d(delta[0] / dir[0], delta[1] / dir[1], delta[2] / dir[2]);
      if (d[0] < d[1] && d[0] < d[2])
      {
        index[0] += int(dir_sign[0]);
      }
      else if (d[1] < d[0] && d[1] < d[2])
      {
        index[1] += int(dir_sign[1]);
      }
      else
      {
        index[2] += int(dir_sign[2]);
      }
    }

    if (progress)
    {
      progress->increment();
    }
  };

#if RAYLIB_PARALLEL_GRID
  tbb::parallel_for<unsigned>(0u, unsigned(cloud.rayCount()), add_ray);
#else   // RAYLIB_PARALLEL_GRID
  const unsigned int count = static_cast<unsigned int>(cloud.rayCount());
//  #pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < count; ++i)
  {
    add_ray(i);
  }
#endif  // RAYLIB_PARALLEL_GRID
}

double Merger::voxelSizeForCloud(const Cloud &cloud) const
{
  double voxel_size = config_.voxel_size;
  if (voxel_size <= 0)
  {
    if (cloud.rayCount() > 0)
    {
      voxel_size = (config_.voxel_size > 0) ? config_.voxel_size : 4.0 * cloud.estimatePointSpacing();
    }
    else
    {
      // Set a reasonable default.
      voxel_size = 0.25;
    }
  }
  return voxel_size;
}

void Merger::markIntersectedEllipsoids(const Cloud &cloud, const Grid<unsigned> &ray_grid,
                                       std::vector<Bool> *transient_ray_marks, double num_rays, bool self_transient,
                                       Progress *progress, bool ellipsoid_cloud_first)
{
  progress->begin("transient-mark-ellipsoids", cloud.rayCount());

  // Check each ellipsoid against the ray grid for intersections.
#if RAYLIB_WITH_TBB
  // Declare thread local for ellipsoid marking
  using ThreadLocalRayMarkers = tbb::enumerable_thread_specific<EllipsoidTransientMarker>;
  ThreadLocalRayMarkers thread_markers(EllipsoidTransientMarker(cloud.rayCount()));

  auto tbb_process_ellipsoid = [this, &cloud, &ray_grid, transient_ray_marks, &num_rays, &thread_markers,
                                ellipsoid_cloud_first, progress, self_transient](size_t ellipsoid_id)  //
  {
    // Resolve the ray marker for this thread.
    EllipsoidTransientMarker &marker = thread_markers.local();
    marker.mark(&ellipsoids_[ellipsoid_id], transient_ray_marks, cloud, ray_grid, num_rays, config_.merge_type,
                self_transient, ellipsoid_cloud_first);
    progress->increment();
  };
  tbb::parallel_for<size_t>(0u, cloud.rayCount(), tbb_process_ellipsoid);
#else   // RAYLIB_WITH_TBB
  std::vector<bool> ray_tested;
  ray_tested.resize(cloud.rayCount(), false);
  EllipsoidTransientMarker ellipsoid_maker(cloud.rayCount());
 // #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < ellipsoids_.size(); ++i)
  {
    ellipsoid_maker.mark(&ellipsoids_[i], transient_ray_marks, cloud, ray_grid, num_rays, config_.merge_type,
                         self_transient, ellipsoid_cloud_first);
    progress->increment();
  }
#endif  // RAYLIB_WITH_TBB
}


void Merger::finaliseFilter(const Cloud &cloud, const std::vector<Bool> &transient_ray_marks)
{
  // Lastly, generate the new ray clouds from this sphere information
  for (size_t i = 0; i < ellipsoids_.size(); i++)
  {
    RGBA col = cloud.colours[i];
    if (config_.colour_cloud)
    {
      col.red = (uint8_t)0;
      col.blue = (uint8_t)(ellipsoids_[i].opacity * 255.0);
      col.green = (uint8_t)((double)ellipsoids_[i].num_gone / ((double)ellipsoids_[i].num_gone + 10.0) * 255.0);
    }

    if (ellipsoids_[i].transient || transient_ray_marks[i])
    {
      difference_.starts.emplace_back(cloud.starts[i]);
      difference_.ends.emplace_back(cloud.ends[i]);
      difference_.times.emplace_back(cloud.times[i]);
      difference_.colours.emplace_back(col);
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
}  // namespace ray
