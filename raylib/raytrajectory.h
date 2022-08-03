// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYTRAJECTORY_H
#define RAYLIB_RAYTRAJECTORY_H

#include "raylib/raylibconfig.h"

#include "raypose.h"
#include "rayutils.h"

namespace ray
{
/// Simple class to represent a moving point over time
class RAYLIB_EXPORT Trajectory
{
public:
  inline std::vector<Eigen::Vector3d> &points() { return points_; }
  inline const std::vector<Eigen::Vector3d> &points() const { return points_; }
  inline std::vector<double> &times() { return times_; }
  inline const std::vector<double> &times() const { return times_; }

  /// Save trajectory to a text file. One line per Node
  bool save(const std::string &file_name);

  /// Load trajectory from file. The file is expected to be a text file, with one Node entry per line
  bool load(const std::string &file_name);

  /// Interpolation of the set @c starts based on the @c times_ of the trajectory
  void calculateStartPoints(const std::vector<double> &times, std::vector<Eigen::Vector3d> &starts);

  /// Nearest position node on the trajectory to the given @c time
  Eigen::Vector3d nearest(double time) const;

  /// Linear interpolation/extrapolation of nearest neighbours at given time.
  /// If 'extrapolate' is false, outlier times will clamp to the start or end value
  Eigen::Vector3d linear(double time, bool extrapolate = true) const;

private:
  inline size_t getIndexAndNormaliseTime(double &time) const
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
  std::vector<Eigen::Vector3d> points_;
  std::vector<double> times_;
};

/// Basic structure to represent a single node in a trajectory
struct RAYLIB_EXPORT TrajectoryNode
{
  Eigen::Vector3d point;
  double time;
};
/// Save a trajectory, using the alternate, node representation
bool RAYLIB_EXPORT saveTrajectory(const std::vector<TrajectoryNode> &nodes, const std::string &file_name);
}  // namespace ray

#endif  // RAYLIB_RAYTRAJECTORY_H
