// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrajectory.h"

namespace ray
{
void Trajectory::calculateStartPoints(const std::vector<double> &times, std::vector<Eigen::Vector3d> &starts)
{
  if (points_.empty() || times_.empty())
    std::cout << "Warning: can only calculate start points when a trajectory is available" << std::endl;

  starts.resize(times.size());
  for (size_t i = 0; i < times.size(); i++) starts[i] = linear(times[i]);
}

bool Trajectory::save(const std::string &file_name)
{
  std::cout << "saving trajectory " << file_name << std::endl;
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (ofs.fail())
  {
    std::cerr << "error: cannot open file " << file_name << " for writing" << std::endl;
    return false;
  }
  ofs.unsetf(std::ios::floatfield);
  ofs.precision(15);
  ofs << "%time x y z userfields" << std::endl;
  for (size_t i = 0; i < points_.size(); i++)
  {
    const Eigen::Vector3d &pos = points_[i];
    ofs << times_[i] << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << std::endl;
  }
  return true;
}

bool saveTrajectory(const std::vector<TrajectoryNode> &nodes, const std::string &file_name)
{
  std::cout << "saving trajectory " << file_name << std::endl;
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (ofs.fail())
  {
    std::cerr << "error: cannot open file " << file_name << " for writing" << std::endl;
    return false;
  }
  ofs.unsetf(std::ios::floatfield);
  ofs.precision(15);
  ofs << "%time x y z userfields" << std::endl;
  for (auto &node : nodes)
  {
    ofs << node.time << " " << node.point[0] << " " << node.point[1] << " " << node.point[2] << " " << std::endl;
  }
  return true;
}


/**Loads the trajectory into the supplied vector and returns if successful*/
bool Trajectory::load(const std::string &file_name)
{
  std::cout << "loading trajectory " << file_name << std::endl;
  std::string line;
  int size = 0;
  {
    std::ifstream ifs(file_name.c_str(), std::ios::in);
    if (!ifs)
    {
      std::cerr << "Failed to open trajectory file: " << file_name << std::endl;
      return false;
    }
    ASSERT(!ifs.fail());

    while (getline(ifs, line))
    {
      if (line.length() == 0 || line[0] == '%')
        continue;
      size++;
    }
  }
  std::ifstream ifs(file_name.c_str(), std::ios::in);
  if (!ifs)
  {
    std::cerr << "Failed to open trajectory file: " << file_name << std::endl;
    return false;
  }
  points_.resize(size);
  times_.resize(size);
  bool ordered = true;
  int i = 0; // index along trajectory
  while (i < size)
  {
    if (ifs.fail())
    {
      std::cerr << "Invalid stream when loading trajectory file: " << file_name << std::endl;
      return false;
    }
    getline(ifs, line);
    if (line.length() == 0 || line[0] == '%')
      continue;

    std::istringstream iss(line);
    iss >> times_[i] >> points_[i][0] >> points_[i][1] >> points_[i][2];
    if (iss.fail())
    {
      std::cerr << "Invalid fields at line " << i << " of " << file_name << std::endl;
      return false;
    }
    if (i > 0 && times_[i] < times_[i - 1])
    {
      ordered = false;
    }
    i++;
  }
  if (ifs.fail())
  {
    std::cerr << "Invalid stream when loading trajectory file: " << file_name << std::endl;
    return false;
  }
  if (!ordered)
  {
    std::cout << "Warning: trajectory times not ordered. Ordering them now." << std::endl;

    struct Temp
    {
      double time;
      size_t index;
    };
    std::vector<Temp> time_list(times_.size());
    for (size_t i = 0; i < time_list.size(); i++)
    {
      time_list[i].time = times_[i];
      time_list[i].index = i;
    }
    sort(time_list.begin(), time_list.end(), [](const Temp &a, const Temp &b) { return a.time < b.time; });

    std::vector<Eigen::Vector3d> new_points(points_.size());
    std::vector<double> new_times(times_.size());
    for (size_t i = 0; i < points_.size(); i++)
    {
      new_points[i] = points_[time_list[i].index];
      new_times[i] = times_[time_list[i].index];
      if (!(i % 100))
        std::cout << "time: " << new_times[i] - new_times[0] << std::endl;
    }
    points_ = std::move(new_points);
    times_ = std::move(new_times);
    std::cout << "finished sorting" << std::endl;
  }

  return true;
}

Eigen::Vector3d Trajectory::nearest(double time) const
{
  ASSERT(!points_.empty());
  if (points_.size() == 1)
    return points_[0];
  const size_t index = getIndexAndNormaliseTime(time);
  return time < 0.5 ? points_[index] : points_[index + 1];
}

/// Linear interpolation/extrapolation of nearest neighbours at given time.
/// If 'extrapolate' is false, outlier times will clamp to the start or end value
Eigen::Vector3d Trajectory::linear(double time, bool extrapolate) const
{
  ASSERT(!points_.empty());
  if (points_.size() == 1)
    return points_[0];
  const size_t index = getIndexAndNormaliseTime(time);
  if (!extrapolate)
  {
    if (time < 0.0)
      return points_.front();
    else if (time > 1.0)
      return points_.back();
  }
  return points_[index] * (1 - time) + points_[index + 1] * time;
}
}  // namespace ray