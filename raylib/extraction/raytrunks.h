// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBRANCHES_H
#define RAYLIB_RAYBRANCHES_H

#include "raylib/raylibconfig.h"
#include "raytrunk.h"
#include "../rayutils.h"
#include "../raycloud.h"
#include <map>
#include <iostream>

namespace ray
{
struct Trunks
{
  Trunks(const Cloud &cloud, double midRadius, bool verbose, bool exclude_passing_rays);
  bool save(const std::string &filename);
  static std::vector<std::pair<Eigen::Vector3d, double> > load(const std::string &filename);
  std::vector<Trunk> best_trunks;
};

struct IntegerVoxels
{
  IntegerVoxels(double width, const Eigen::Vector3d offset) : voxel_width(width), offset(offset) {}

  inline Eigen::Vector3i getIndex(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d ind = (pos - offset) / voxel_width;
    return Eigen::Vector3i(int(std::floor(ind[0])), int(std::floor(ind[1])), int(std::floor(ind[2])));
  }
  inline void increment(const Eigen::Vector3d &pos)
  {
    increment(getIndex(pos));
  }
  inline void increment(const Eigen::Vector3i &index)
  {
    auto it = voxel_map.find(index);
    if (it == voxel_map.end())
    {
      voxel_map[index] = 1;
    }
    else
    {
      it->second++;
    }
  }
  inline int get(const Eigen::Vector3i &index) const
  {
    auto it = voxel_map.find(index);
    if (it == voxel_map.end())
      return 0;
    else
      return it->second;
  }
  void forEach(std::function<void(const struct IntegerVoxels &voxels, double width, const Eigen::Vector3d &offset, const Eigen::Vector3i &index, int count)> func)
  {
    for (auto &voxel: voxel_map)
    {
      func(*this, voxel_width, offset, voxel.first, voxel.second);
    }
  }

  std::map<Eigen::Vector3i, int, Vector3iLess> voxel_map;
  double voxel_width;
  Eigen::Vector3d offset;
};

template<typename T>
void writePlainOldData(std::ofstream &out, const T &t)
{
  out.write(reinterpret_cast<const char*>(&t), sizeof(T));
}

template<typename T>
void readPlainOldData(std::ifstream &in, T &t)
{
  in.read(reinterpret_cast<char*>(&t), sizeof(T));
}

template<typename T>
void writePlainOldDataArray(std::ofstream &out, const std::vector<T> &array)
{
  unsigned int size = (unsigned int) array.size();
  out.write(reinterpret_cast<char*>(&size), sizeof(unsigned int));
  for (unsigned int i = 0; i<size; i++)
    writePlainOldData(out, array[i]);
}

template<typename T>
void readPlainOldDataArray(std::ifstream &in, std::vector<T> &array)
{
  unsigned int size;
  in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
  array.resize(size);
  for (unsigned int i = 0; i<size; i++)
    readPlainOldData(in, array[i]);
}


} // namespace ray
#endif // RAYLIB_RAYBRANCHES_H