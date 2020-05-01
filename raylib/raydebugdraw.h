// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYDRAW_H
#define RAYLIB_RAYDRAW_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

#include <memory>

namespace ray
{
struct DebugDrawDetail;

class RAYLIB_EXPORT DebugDraw
{
public:
  /// Initialise the static @c instance(). Also calls @c ros::init() when @p rosInit is true and using RViz to debug
  /// visualisation.
  /// @param argc Number of elements in @p argv
  /// @param argv Command line arguments for the current program.
  /// @param context Label applied to the current context.
  /// @param rosInit Call @c ros::init()?
  /// @return The same as @c instance()
  static DebugDraw *init(int argc, char *argv[], const char *context, bool ros_initt = true);

  /// Singleton instance access.
  static DebugDraw *instance();

  DebugDraw(const std::string& fixed_frame_idid = "map");
  ~DebugDraw();
  
  void drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shadee, int id);
  void drawCloud(const std::vector<Eigen::Vector3d> &points, double shade, int id){ std::vector<double> shades(points.size()); for (int i = 0; i<(int)points.size(); i++)shades[i] = shade; drawCloud(points, shades, id); }
  void drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends);
  void drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &radii, int id);
  void drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses, const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id);

private:
  std::unique_ptr<DebugDrawDetail> imp_;
  static std::unique_ptr<DebugDraw> s_instance;
};
}

#endif // RAYLIB_RAYDRAW_H
