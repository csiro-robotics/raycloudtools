// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTRAJECTORY_H
#define RAYLIB_RAYTRAJECTORY_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"

namespace ray
{

/// Simple class to represent a moving point over time
class RAYLIB_EXPORT Trajectory
{
public:
  std::vector<Eigen::Vector3d> points_;
  std::vector<double> times_;

  /// Save trajectory to a text file. One line per Node
  void save(const std::string &file_name);

  /// Load trajectory from file. The file is expected to be a text file, with one Node entry per line
  bool load(const std::string &file_name);

  void calculateStartPoints(const std::vector<double> &times, std::vector<Eigen::Vector3d> &starts);


  /// Nearest value at a given time
  Eigen::Vector3d nearest(double time) const
  {
    ASSERT(!points_.empty());
    if (points_.size() == 1)
      return points_[0];
    size_t index = getIndexAndNormaliseTime(time);
    return time < 0.5 ? points_[index] : points_[index + 1];
  }
  
  /// Linear interpolation/extrapolation of nearest neighbours at given time.
  /// If 'extrapolate' is false, outlier times will clamp to the start or end value
  Eigen::Vector3d linear(double time, bool extrapolate = true) const
  {
    ASSERT(!points_.empty());
    if (points_.size() == 1)
      return points_[0];
    size_t index = getIndexAndNormaliseTime(time);
    if (!extrapolate)
    {
      if (time < 0.0)
        return points_.front();
      else if (time > 1.0)
        return points_.back();
    }
    return points_[index] * (1-time) + points_[index+1] * time;
  }

protected:
  size_t getIndexAndNormaliseTime(double &time) const
  {
    size_t index = std::lower_bound(times_.begin(), times_.end(), time) - times_.begin();
    if (index == 0)
      index++;
    if (index == times_.size())
      index--;
    index--;
    time = (time - times_[index]) / (times_[index + 1] - times_[index]);
    return index;
  }
};

}

#endif // RAYLIB_RAYTRAJECTORY_H
