// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
//#define HASH_LOOKUP

namespace RAY
{
#if defined HASH_LOOKUP
template<class T> 
struct Grid
{
  Grid(){}
  Grid(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth){ init(boxMin, boxMax, voxelWidth); }
  void init(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth)
  {
    this->boxMin = boxMin;
    this->boxMax = boxMax;
    this->voxelWidth = voxelWidth;
    Eigen::Vector3d diff = (boxMax - boxMin)/voxelWidth;
    dims = Eigen::Vector3i(ceil(diff[0]), ceil(diff[1]), ceil(diff[2]));   

    int bucketSize = dims[0]*dims[1]; // let's assume a surface for the map. This is still significantly better on memory than a 3D grid (dims[0]*dims[1]*dims[2])
    buckets.resize(bucketSize);
    nullCell.index = Eigen::Vector3i(-1,-1,-1);
  }
  struct Cell
  {
    std::vector<T> data;
    Eigen::Vector3i index;
  };
  
  Cell &cell(int x, int y, int z)
  {
    Eigen::Vector3i index(x,y,z);
    int hash = (x*17 + y*101 + z*797) % buckets.size();
    std::vector<Cell> &cells = buckets[hash].cells;
    for (auto &c: cells)
      if (c.index == index)
        return c;
    return nullCell;    
  }
  void insert(int x, int y, int z, const T &value)
  {
    Eigen::Vector3i index(x,y,z);
    int hash = (x*17 + y*101 + z*797) % buckets.size();
    std::vector<Cell> &cellList = buckets[hash].cells;
    for (auto &c: cellList)
    {
      if (c.index == index)
      {
        c.data.push_back(value);
        return;
      }
    }
    Cell newCell;
    newCell.index = index;
    cellList.push_back(newCell);
    cellList.back().data.push_back(value);    
  }

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
  Eigen::Vector3i dims;

  void report()
  {
    int count = 0;
    int totalCount = 0;
    int dataCount = 0;
    for (auto &bucket: buckets)
    {
      int size= bucket.cells.size();
      if (size>0)
        count++;
      totalCount += size;
      for (auto &cell: bucket.cells)
        dataCount += cell.data.size();
    }
    std::cout << "buckets filled: " << count << " out of " << buckets.size() << " buckets, which is " << 100.0*(double)count/(double)buckets.size() << "\%" << std::endl;
    std::cout << "voxels filled: " << totalCount << ", average size of filled bucket: " << (double)totalCount/(double)count << std::endl;
    std::cout << "average data per filled voxel: " << (double)dataCount/(double)totalCount << std::endl;
    std::cout << "total data stored: " << dataCount << std::endl;
  }

protected:
  struct Bucket
  {
    std::vector<Cell> cells;
  };
  std::vector<Bucket> buckets;
  Cell nullCell;
};

#else

template<class T> 
struct Grid
{
  Grid(){}
  Grid(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth){ init(boxMin, boxMax, voxelWidth); }
  void init(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth)
  {
    this->boxMin = boxMin;
    this->boxMax = boxMax;
    this->voxelWidth = voxelWidth;
    Eigen::Vector3d diff = (boxMax - boxMin)/voxelWidth;
    dims = Eigen::Vector3i(ceil(diff[0]), ceil(diff[1]), ceil(diff[2]));   
    cells.resize(dims[0]*dims[1]*dims[2]);
  }
  struct Cell
  {
    std::vector<T> data;
  };
  inline Cell &cell(int x, int y, int z)
  {
    if (x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2])
      return nullCell;
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }
  inline void insert(int x, int y, int z, const T &value)
  {
    if (x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2])
      std::cout << "warning: bad input coordinates: " << x << ", " << y << ", " << z << std::endl;
    cell(x, y, z).data.push_back(value);
  }
  void report()
  {
    int count = 0;
    int totalCount = 0;

    for (auto &cell: cells)
    {
      int size = cell.data.size();
      if (size > 0)
        count ++;
      totalCount += size;
    }
    std::cout << "voxels filled: " << count << " out of " << cells.size() << " cells, which is " << 100.0*(double)count/(double)cells.size() << "\%" << std::endl;
    std::cout << "average data per filled voxel: " << (double)totalCount/(double)count << std::endl;
    std::cout << "total data stored: " << totalCount << std::endl;
  }

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
  Eigen::Vector3i dims;
protected:
  std::vector<Cell> cells;
  Cell nullCell;
};
#endif

}
