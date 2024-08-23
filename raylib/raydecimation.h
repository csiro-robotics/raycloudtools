// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYDECIMATION_H
#define RAYLIB_RAYDECIMATION_H

#include <iostream>
#include <limits>
#include "raycloud.h"
#include "rayutils.h"

namespace ray
{
/// @brief subsample to 1 point per @c vox_width wide voxel in metres
/// This is a spatially even subsampling, but also emphasises outlier as a side-effect
bool RAYLIB_EXPORT decimateSpatial(const std::string &file_stub, double vox_width);

/// @brief subsample to every @c num_rays rays
/// This is an unbiased subsampling, but will be over-sampled in stationary areas as a side-effect
/// Note that while this is called temporal decimation, it decimates evenly in file order, which isn't 
/// necessarily temporal order. Though it typically is stored that way on single scans.
bool RAYLIB_EXPORT decimateTemporal(const std::string &file_stub, int num_rays);

/// @brief subsample to @c num_rays rays (temporally decimated) for each @c vox_width wide voxel
/// This allows a more even distribution of points while maintaining details better than pure spatial decimation
bool RAYLIB_EXPORT decimateSpatioTemporal(const std::string &file_stub, double vox_width, int num_rays);

/// @brief Maintains a maximum number of rays intersecting each voxel. This has some ambiguity, but is a useful routine
/// as it maintains the integrity of the full ray cloud including free space, so is better for combine operations
/// By contrast, standard spatial decimation removes free space whenever the end points coincide 
bool RAYLIB_EXPORT decimateRaysSpatial(const std::string &file_stub, double vox_width);


/// @brief decimate to no more than 1 point per voxel of width @c radius_per_length x ray length. 
/// This is used when error is proportional to ray length, prioritising closer measurements and leaving distant areas sparse
bool RAYLIB_EXPORT decimateAngular(const std::string &file_stub, double radius_per_length);


struct Subsampler
{
  inline bool operator()(const Eigen::Vector3i &p, int i)
  {
    if (voxel_set.insert(p).second)
    {
      subsample.push_back(i);
      return true;
    }
    return false;
  }
  std::vector<int64_t> subsample;
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
};

inline int sign(double x)
{
  return (x > 0.0) - (x < 0.0);
}

// for similar appraoch see: https://github.com/StrandedKitty/tiles-intersect/blob/master/src/index.js
template<class T> 
void walkGrid(const Eigen::Vector3d &start, const Eigen::Vector3d &end, T &object, int i)
{
  Eigen::Vector3i p = Eigen::Vector3d(std::floor(start[0]), std::floor(start[1]), std::floor(start[2])).cast<int>();
  
  if (object(p, i))
  {
    return; // only adding to one cell
  }
  
  Eigen::Vector3i step(sign(end[0] - start[0]), sign(end[1] - start[1]), sign(end[2] - start[2]));
  Eigen::Vector3d tmax, tdelta;
  for (int j = 0; j<3; j++)
  {
    step[j] = sign(end[j] - start[j]);
    double to = std::abs(start[j] - p[j] - (double)std::max(0, step[j]));        
    double dir = std::max(std::numeric_limits<double>::min(), std::abs(start[j] - end[j]));
    tmax[j] = to / dir;
    tdelta[j] = 1.0 / dir;
  }
  
  Eigen::Vector3i target = Eigen::Vector3d(std::floor(end[0]), std::floor(end[1]), std::floor(end[2])).cast<int>();
  while (p != target) 
  {
    int ax = tmax[0] < tmax[1] && tmax[0] < tmax[2] ? 0 : (tmax[1] < tmax[2] ? 1 : 2);
    p[ax] += step[ax];
    tmax[ax] += tdelta[ax];
    if (object(p, i))
    {
      break; // only adding to one cell
    }          
  }     
}
}  // namespace ray

#endif  // RAYLIB_RAYDECIMATION_H
