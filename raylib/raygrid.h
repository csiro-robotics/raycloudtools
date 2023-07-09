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

#if RAYLIB_WITH_TBB
#define RAYLIB_PARALLEL_GRID 1
#if RAYLIB_PARALLEL_GRID
#include <tbb/spin_mutex.h>
#endif  // RAYLIB_PARALLEL_GRID
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

/// 3D grid container class based on hash lookup, to accelerate the access to spatial data by location
/// A hash lookup is used because ray cloud geometry is generally sparse, and so continuous 3D voxel arrays are memory
/// intensive
template <class T>
class Grid
{
public:
#if RAYLIB_PARALLEL_GRID
  using Mutex = tbb::spin_mutex;
#endif  // RAYLIB_PARALLEL_GRID

  class Cell
  {
  public:
    std::vector<T> data;
    Eigen::Vector3i index;
#if RAYLIB_PARALLEL_GRID
    Mutex mutex;
#endif  // RAYLIB_PARALLEL_GRID

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

      coord(i) = (coord(i) - box_min(i)) / voxel_width;
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

  /// the grid is axis aligned, so initialised from a bounding box and a voxel width
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
    // linearly search through the bucket
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
    insert(index, value);
  }

  void addCell(const Eigen::Vector3i &index)
  {
    int hash = hashFunc(index[0], index[1], index[2]);
    Bucket &bucket = buckets_.at(hash);
#if RAYLIB_PARALLEL_GRID
    Mutex::scoped_lock bucket_lock(bucket.mutex);
#endif  // RAYLIB_PARALLEL_GRID
    for (auto &c : bucket.cells)
    {
      if (c.index == index)
      {
        return;
      }
    }
    Cell new_cell;
    new_cell.index = index;
    bucket.cells.emplace_back(new_cell); // fill with empty cell
  }

  void insert(const Eigen::Vector3i &index, const T &value)
  {
    int hash = hashFunc(index[0], index[1], index[2]);
    Bucket &bucket = buckets_.at(hash);
#if RAYLIB_PARALLEL_GRID
    Mutex::scoped_lock bucket_lock(bucket.mutex);
#endif  // RAYLIB_PARALLEL_GRID
    for (auto &c : bucket.cells)
    {
      if (c.index == index)
      {
#if RAYLIB_PARALLEL_GRID
        Mutex::scoped_lock cell_lock(c.mutex);
#endif  // RAYLIB_PARALLEL_GRID
        if (c.index == index)
        {
          c.data.emplace_back(value);
          return;
        }
      }
    }
    bucket.cells.emplace_back(Cell(index, value));
  }

  // only inserts into a cell that exists
  void insertIfCellExists(const Eigen::Vector3i &index, const T &value)
  {
    int hash = hashFunc(index[0], index[1], index[2]);
    Bucket &bucket = buckets_.at(hash);
#if RAYLIB_PARALLEL_GRID
    Mutex::scoped_lock bucket_lock(bucket.mutex);
#endif  // RAYLIB_PARALLEL_GRID
    for (auto &c : bucket.cells)
    {
      if (c.index == index)
      {
#if RAYLIB_PARALLEL_GRID
        Mutex::scoped_lock cell_lock(c.mutex);
#endif  // RAYLIB_PARALLEL_GRID
        if (c.index == index)
        {
          c.data.emplace_back(value);
          return;
        }
      }
    }
  }

  /// debugging statistics on the grid structure. This can be used to assess how efficient this grid
  /// structure is for a given @c voxel_width.
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

  /// applies the @c visit function for all cells in the grid
  void walkCells(const WalkCellsVisitFunction &visit) const
  {
    for (const auto &bucket : buckets_)
    {
      for (const auto &cell : bucket.cells)
      {
        visit(*this, cell);
      }
    }
  }

  Eigen::Vector3d box_min, box_max;
  double voxel_width;
  Eigen::Vector3i dims;

protected:
  class Bucket
  {
  public:
    std::vector<Cell> cells;
#if RAYLIB_PARALLEL_GRID
    Mutex mutex;
#endif  // RAYLIB_PARALLEL_GRID

    inline Bucket() = default;
    inline Bucket(const Bucket &other)
      : cells(other.cells)
    {}
    inline Bucket(Bucket &&other)
      : cells(std::move(other.cells))
    {}
  };

  inline int hashFunc(int x, int y, int z) const { return int((x * 17 + y * 101 + z * 797) % buckets_.size()); }

  std::vector<Bucket> buckets_;
  Cell null_cell_;
};

template <class T>
class ContiguousGrid
{
public:
  ContiguousGrid() {}
  ContiguousGrid(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width,
                 double voxel_height = 0)
  {
    init(box_min, box_max, voxel_width, voxel_height);
  }
  void init(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width, double voxel_height = 0)
  {
    this->box_min = box_min;
    this->voxel_width = voxel_width;
    this->voxel_height = voxel_height ? voxel_height : voxel_width;
    Eigen::Vector3d diff = (box_max - box_min) / voxel_width;
    dims = Eigen::Vector3i(ceil(diff[0]), ceil(diff[1]), ceil((box_max[2] - box_min[2]) / voxel_height));
    cells.resize(dims[0] * dims[1] * dims[2]);
  }
  struct Cell
  {
    std::vector<T> data;
  };
  Eigen::Vector3i index(const Eigen::Vector3d &spatial_pos) const
  {
    Eigen::Vector3d coord = spatial_pos - box_min;
    coord[0] /= voxel_width;
    coord[1] /= voxel_width;
    coord[2] /= voxel_height;
    return coord.cast<int>();
  }
  inline Cell &cell(const Eigen::Vector3i &index)
  {
    const double &x = index[0];
    const double &y = index[1];
    const double &z = index[2];
    if (x < 0 || x >= dims[0] || y < 0 || y >= dims[1] || z < 0 || z >= dims[2])
      return null_cell_;
    return cells[x + dims[0] * y + dims[0] * dims[1] * z];
  }
  inline const Cell &cell(const Eigen::Vector3i &index) const
  {
    const double &x = index[0];
    const double &y = index[1];
    const double &z = index[2];
    if (x < 0 || x >= dims[0] || y < 0 || y >= dims[1] || z < 0 || z >= dims[2])
      return null_cell_;
    return cells[x + dims[0] * y + dims[0] * dims[1] * z];
  }
  inline void insert(const Eigen::Vector3i &index, const T &value)
  {
    const double &x = index[0];
    const double &y = index[1];
    const double &z = index[2];
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
              << 100.0 * (double)count / (double)cells.size() << "%" << std::endl;
    std::cout << "average data per filled voxel: " << (double)totalCount / (double)count << std::endl;
    std::cout << "total data stored: " << totalCount << std::endl;
  }

  Eigen::Vector3d box_min;
  double voxel_width, voxel_height;
  Eigen::Vector3i dims;

protected:
  std::vector<Cell> cells;
  Cell null_cell_;
};

}  // namespace ray

#endif  // RAYLIB_RAYGRID_H
