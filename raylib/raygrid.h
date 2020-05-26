// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe, Kazys Stepanas
#ifndef RAYLIB_RAYGRID_H
#define RAYLIB_RAYGRID_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

#include <functional>

#define HASH_LOOKUP

#if RAYLIB_WITH_TBB || defined(HASH_LOOKUP)
#define RALIB_PARALLEL_GRID 1
#if RALIB_PARALLEL_GRID
#include <tbb/spin_mutex.h>
#endif  // RALIB_PARALLEL_GRID
#endif  // RAYLIB_WITH_TBB

namespace ray
{
struct RAYLIB_EXPORT GridRayInfo
{
  Eigen::Vector3d ray_start;
  Eigen::Vector3d ray_end;
  Eigen::Vector3d ray_direction;
  double ray_length;
};

#if defined HASH_LOOKUP
template <class T>
struct Grid
{
#if RALIB_PARALLEL_GRID
  using Mutex = tbb::spin_mutex;
#endif  // RALIB_PARALLEL_GRID

  struct Cell
  {
    std::vector<T> data;
    Eigen::Vector3i index;
#if RALIB_PARALLEL_GRID
    Mutex mutex;
#endif  // RALIB_PARALLEL_GRID

    inline Cell() {}
    inline Cell(const Eigen::Vector3i &index, T initial_datum)
      : index(index)
    {
      data.push_back(initial_datum);
    }

    /// Copy constructor, excludes mutex
    inline Cell(const Cell &other)
      : data(other.data)
      , index(other.index)
    {}

    /// RValue constructor, excludes mutex
    inline Cell(Cell &&other)
      : data(std::move(other.data))
      , index(other.index)
    {}
  };

  using WalkVoxelsVisitFunction =
    std::function<void(const Grid<T> &, const Eigen::Vector3i &, const GridRayInfo &info)>;
  using WalkCellsVisitFunction = std::function<void(const Grid<T> &, const Cell &)>;

  Grid() {}
  Grid(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width)
  {
    init(box_min, box_max, voxel_width);
  }

  /// Generate a grid indexer from a spatial position.
  /// @p spatial_pos is first clamped to the bounds if @c clamp is true. Otherwise, the return value may be out of
  /// range.
  Eigen::Vector3i index(const Eigen::Vector3d &spatial_pos, bool clamp = true) const
  {
    Eigen::Vector3d coord = spatial_pos;
    // Bounds check.
    for (int i = 0; i < 3; ++i)
    {
      if (coord(i) < box_min(i))
      {
        if (!clamp)
        {
          return Eigen::Vector3i(-1, -1, -1);
        }
        coord(i) = box_min(i);
      }
      if (coord(i) > box_max(i))
      {
        if (!clamp)
        {
          return Eigen::Vector3i(-1, -1, -1);
        }
        coord(i) = box_max(i);
      }

      coord(i) = (coord(i) - box_min(i)) / voxel_width + 0.5;
    }

    return coord.cast<int>();
  }

  Eigen::Vector3d voxelCentre(const Eigen::Vector3i &index) const
  {
    return voxelCentre(index.x(), index.y(), index.z());
  }

  Eigen::Vector3d voxelCentre(int x, int y, int z) const
  {
    return Eigen::Vector3d(box_min.x() + x * voxel_width + 0.5 * voxel_width,
                           box_min.y() + y * voxel_width + 0.5 * voxel_width,
                           box_min.z() + z * voxel_width + 0.5 * voxel_width);
  }

  void init(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width)
  {
    this->box_min = box_min;
    this->box_max = box_max;
    this->voxel_width = voxel_width;
    Eigen::Vector3d diff = (box_max - box_min) / voxel_width;
    dims = Eigen::Vector3i(diff.array().ceil().cast<int>());

    int bucket_size = dims[0] * dims[1];  // let's assume a surface for the map. This is still significantly better on
                                          // memory than a 3D grid (dims[0]*dims[1]*dims[2])
    buckets_.resize(bucket_size);
    null_cell_.index = Eigen::Vector3i(-1, -1, -1);
  }

  Cell &cell(int x, int y, int z) { return cell(Eigen::Vector3i(x, y, z)); }
  Cell &cell(const Eigen::Vector3i &index)
  {
    int hash = hashFunc(index.x(), index.y(), index.z());
    Bucket &bucket = buckets_.at(hash);
    for (auto &c : bucket.cells)
      if (c.index == index)
        return c;
    return null_cell_;
  }


  const Cell &cell(int x, int y, int z) const { return cell(Eigen::Vector3i(x, y, z)); }
  const Cell &cell(const Eigen::Vector3i &index) const
  {
    int hash = hashFunc(index.x(), index.y(), index.z());
    const Bucket &bucket = buckets_.at(hash);
    for (const auto &c : bucket.cells)
      if (c.index == index)
        return c;
    return null_cell_;
  }

  void insert(int x, int y, int z, const T &value)
  {
    Eigen::Vector3i index(x, y, z);
    int hash = hashFunc(x, y, z);
    Bucket &bucket = buckets_.at(hash);
#if RALIB_PARALLEL_GRID
    Mutex::scoped_lock bucket_lock(bucket.mutex);
#endif  // RALIB_PARALLEL_GRID
    for (auto &c : bucket.cells)
    {
      if (c.index == index)
      {
#if RALIB_PARALLEL_GRID
        Mutex::scoped_lock cell_lock(c.mutex);
#endif  // RALIB_PARALLEL_GRID
        if (c.index == index)
        {
          c.data.emplace_back(value);
          return;
        }
      }
    }
    bucket.cells.emplace_back(Cell(index, value));
  }

  void report()
  {
    size_t count = 0;
    size_t total_count = 0;
    size_t data_count = 0;
    for (auto &bucket : buckets_)
    {
      size_t size = bucket.cells.size();
      if (size > 0)
        count++;
      total_count += size;
      for (auto &cell : bucket.cells)
      {
        data_count += cell.data.size();
      }
    }
    std::cout << "buckets filled: " << count << " out of " << buckets_.size() << " buckets, which is "
              << 100.0 * (double)count / (double)buckets_.size() << "%%" << std::endl;
    std::cout << "voxels filled: " << total_count
              << ", average size of filled bucket: " << (double)total_count / (double)count << std::endl;
    std::cout << "average data per filled voxel: " << (double)data_count / (double)total_count << std::endl;
    std::cout << "total data stored: " << data_count << std::endl;
  }

  size_t walkVoxels(const Eigen::Vector3d &start, const Eigen::Vector3d &end, const WalkVoxelsVisitFunction &visit,
                    bool include_end_point, double length_epsilon = 1e-6) const;

  void walkCells(const WalkCellsVisitFunction &visit) const
  {
    for (const auto &bucket : buckets_)
    {
      for (const auto &cell : bucket.cells)
      {
        visit(*this, cell);
      }
    }
    visit(*this, null_cell_);
  }

  Eigen::Vector3d box_min, box_max;
  double voxel_width;
  Eigen::Vector3i dims;

protected:
  struct Bucket
  {
    std::vector<Cell> cells;
#if RALIB_PARALLEL_GRID
    Mutex mutex;
#endif  // RALIB_PARALLEL_GRID

    inline Bucket() = default;
    inline Bucket(const Bucket &other) : cells(other.cells) {}
    inline Bucket(Bucket &&other) : cells(std::move(other.cells)) {}
  };

  inline int hashFunc(int x, int y, int z) const { return int((x * 17 + y * 101 + z * 797) % buckets_.size()); }

  std::vector<Bucket> buckets_;
  Cell null_cell_;
};


template <typename T>
size_t Grid<T>::walkVoxels(const Eigen::Vector3d &start_point, const Eigen::Vector3d &end_point,
                           const WalkVoxelsVisitFunction &visit, bool include_end_point, double length_epsilon) const
{
  // see "A Faster Voxel Traversal Algorithm for Ray Tracing" by Amanatides & Woo
  // Convert the start and end points to voxel indices. Do not clamp. Clamping would change the angle of the line.
  // FIXME: (KS) support clamping the ray to the bounds.
  Eigen::Vector3i start_index = index(start_point, false);
  Eigen::Vector3i end_index = index(end_point, false);
  Eigen::Vector3d direction = end_point - start_point;
  double length = direction.squaredNorm();

  // Very small segments which straddle a voxel boundary can be problematic. We want to avoid
  // a sqrt on a very small number, but be robust enough to handle the situation.
  // To that end, we skip normalising the direction for directions below a tenth of the voxel.
  // Then we will exit either with start/end voxels being the same, or we will step from start
  // to end in one go.
  const bool valid_length = (length >= length_epsilon * length_epsilon);
  if (valid_length)
  {
    length = std::sqrt(length);
    direction *= 1.0 / length;
  }

  GridRayInfo ray_info;
  ray_info.ray_start = start_point;
  ray_info.ray_end = end_point;
  ray_info.ray_direction = direction;
  ray_info.ray_length = length;

  if (start_index == end_index)
  {
    if (include_end_point)
    {
      visit(*this, end_index, ray_info);
    }
    return 1u;
  }

  if (!valid_length)
  {
    // Start/end points are in different, but adjacent voxels. Prevent issues with the loop by
    // early out.
    visit(*this, start_index, ray_info);
    if (include_end_point)
    {
      visit(*this, end_index, ray_info);
      return 2;
    }
    return 1;
  }

  int step[3] = { 0 };
  Eigen::Vector3d voxel;
  double time_max[3];
  double time_delta[3];
  double time_limit[3];
  double next_voxel_border;
  double direction_axis_inv;
  size_t added = 0;
  Eigen::Vector3i current_index = start_index;

  voxel = voxelCentre(current_index);

  // Compute step direction, increments and maximums along each axis.
  for (unsigned i = 0; i < 3; ++i)
  {
    if (direction[i] != 0)
    {
      direction_axis_inv = 1.0 / direction[i];
      step[i] = (direction[i] > 0) ? 1 : -1;
      // Time delta is the ray time between voxel boundaries calculated for each axis.
      time_delta[i] = voxel_width * std::abs(direction_axis_inv);
      // Calculate the distance from the origin to the nearest voxel edge for this axis.
      next_voxel_border = voxel[i] + step[i] * 0.5 * voxel_width;
      time_max[i] = (next_voxel_border - start_point[i]) * direction_axis_inv;
      time_limit[i] = std::abs((end_point[i] - start_point[i]) * direction_axis_inv);  // +0.5f * voxel_width;
    }
    else
    {
      time_max[i] = time_delta[i] = std::numeric_limits<double>::max();
      time_limit[i] = 0;
    }
  }

  int axis = 0;
  bool limit_reached = false;
  while (!limit_reached && current_index != end_index)
  {
    // Only visit indices which are in range.
    if (0 <= current_index.x() && current_index.x() < dims.x() &&  //
        0 <= current_index.y() && current_index.y() < dims.y() &&  //
        0 <= current_index.z() && current_index.z() < dims.z())
    {
      visit(*this, current_index, ray_info);
    }
    ++added;
    axis = (time_max[0] < time_max[2]) ? ((time_max[0] < time_max[1]) ? 0 : 1) : ((time_max[1] < time_max[2]) ? 1 : 2);
    limit_reached = std::abs(time_max[axis]) > time_limit[axis];
    current_index[axis] += step[axis];
    time_max[axis] += time_delta[axis];
  }

  if (include_end_point)
  {
    visit(*this, end_index, ray_info);
    ++added;
  }

  // assert(added);
  return added;
}

#else   // HASH_LOOKUP

template <class T>
struct Grid
{
  Grid() {}
  Grid(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width)
  {
    init(box_min, box_max, voxel_width);
  }
  void init(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width)
  {
    this->box_min = box_min;
    this->box_max = box_max;
    this->voxel_width = voxel_width;
    Eigen::Vector3d diff = (box_max - box_min) / voxel_width;
    dims = Eigen::Vector3i(ceil(diff[0]), ceil(diff[1]), ceil(diff[2]));
    cells.resize(dims[0] * dims[1] * dims[2]);
  }
  struct Cell
  {
    std::vector<T> data;
  };
  inline Cell &cell(int x, int y, int z)
  {
    if (x < 0 || x >= dims[0] || y < 0 || y >= dims[1] || z < 0 || z >= dims[2])
      return null_cell_;
    return cells[x + dims[0] * y + dims[0] * dims[1] * z];
  }
  inline void insert(int x, int y, int z, const T &value)
  {
    if (x < 0 || x >= dims[0] || y < 0 || y >= dims[1] || z < 0 || z >= dims[2])
      std::cout << "warning: bad input coordinates: " << x << ", " << y << ", " << z << std::endl;
    cell(x, y, z).data.emplace_back(value);
  }
  void report()
  {
    int count = 0;
    int totalCount = 0;

    for (auto &cell : cells)
    {
      int size = cell.data.size();
      if (size > 0)
        count++;
      totalCount += size;
    }
    std::cout << "voxels filled: " << count << " out of " << cells.size() << " cells, which is "
              << 100.0 * (double)count / (double)cells.size() << "\%" << std::endl;
    std::cout << "average data per filled voxel: " << (double)totalCount / (double)count << std::endl;
    std::cout << "total data stored: " << totalCount << std::endl;
  }

  Eigen::Vector3d box_min, box_max;
  double voxel_width;
  Eigen::Vector3i dims;

protected:
  std::vector<Cell> cells;
  Cell null_cell_;
};
#endif  // HASH_LOOKUP

}  // namespace ray

#endif  // RAYLIB_RAYGRID_H
