// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYGRID_H
#define RAYLIB_RAYGRID_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#define HASH_LOOKUP

namespace ray
{
#if defined HASH_LOOKUP
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
    dims = Eigen::Vector3i(diff.array().ceil().cast<int>());

    int bucket_size = dims[0] * dims[1];  // let's assume a surface for the map. This is still significantly better on
                                          // memory than a 3D grid (dims[0]*dims[1]*dims[2])
    buckets_.resize(bucket_size);
    null_cell_.index = Eigen::Vector3i(-1, -1, -1);
  }
  struct Cell
  {
    std::vector<T> data;
    Eigen::Vector3i index;
  };

  Cell &cell(int x, int y, int z)
  {
    Eigen::Vector3i index(x, y, z);
    int hash = (x * 17 + y * 101 + z * 797) % (int)buckets_.size();
    std::vector<Cell> &cells = buckets_[hash].cells;
    for (auto &c : cells)
      if (c.index == index)
        return c;
    return null_cell_;
  }
  void insert(int x, int y, int z, const T &value)
  {
    Eigen::Vector3i index(x, y, z);
    int hash = (x * 17 + y * 101 + z * 797) % (int)buckets_.size();
    std::vector<Cell> &cell_list = buckets_[hash].cells;
    for (auto &c : cell_list)
    {
      if (c.index == index)
      {
        c.data.push_back(value);
        return;
      }
    }
    Cell new_cell;
    new_cell.index = index;
    cell_list.push_back(new_cell);
    cell_list.back().data.push_back(value);
  }

  Eigen::Vector3d box_min, box_max;
  double voxel_width;
  Eigen::Vector3i dims;

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
      for (auto &cell : bucket.cells) data_count += cell.data.size();
    }
    std::cout << "buckets filled: " << count << " out of " << buckets_.size() << " buckets, which is "
              << 100.0 * (double)count / (double)buckets_.size() << "%%" << std::endl;
    std::cout << "voxels filled: " << total_count
              << ", average size of filled bucket: " << (double)total_count / (double)count << std::endl;
    std::cout << "average data per filled voxel: " << (double)data_count / (double)total_count << std::endl;
    std::cout << "total data stored: " << data_count << std::endl;
  }

protected:
  struct Bucket
  {
    std::vector<Cell> cells;
  };
  std::vector<Bucket> buckets_;
  Cell null_cell_;
};

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
    cell(x, y, z).data.push_back(value);
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
