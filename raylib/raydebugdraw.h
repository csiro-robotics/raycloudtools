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

  DebugDraw(const std::string &fixed_frame_idid = "map");
  ~DebugDraw();

  /// Draw point cloud, each defined by a position and a shade. The id defines the index
  void drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &point_shade, int id);

  /// Draw point cloud, defined by a set of points and a single shade and id
  void drawCloud(const std::vector<Eigen::Vector3d> &points, double shade, int id)
  {
    std::vector<double> shades(points.size());
    for (int i = 0; i < (int)points.size(); i++) shades[i] = shade;
    drawCloud(points, shades, id);
  }

  // Draw lines, as defined by their start, end and colour.
  void drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                 const std::vector<Eigen::Vector3d> &colours = std::vector<Eigen::Vector3d>());

  // Draw list of cylinders, defined by their start, end, radius and colour. The id allows multiple sets to be drawn independently, where each new call
  // to drawCylinders only replaces the cylinders of the chosen id. 
  void drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                     const std::vector<double> &radii, int id, const std::vector<Eigen::Vector4d> &colours = std::vector<Eigen::Vector4d>());
  
  // Draw ellipsoids, defined by their centre, pose, colour and vector of radii. The ID allows independent sets
  void drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses,
                      const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id);

private:
  std::unique_ptr<DebugDrawDetail> imp_;
  static std::unique_ptr<DebugDraw> s_instance;
};
}  // namespace ray

#endif  // RAYLIB_RAYDRAW_H
