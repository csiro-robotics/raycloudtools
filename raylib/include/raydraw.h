// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#if defined USE_ROS
#include <ros/ros.h>
#endif

namespace RAY
{
struct DebugDraw
{
  DebugDraw(const std::string& fixedFrameId = "map");
  
  void drawCloud(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &pointShade, int id);
  void drawCloud(const std::vector<Eigen::Vector3d> &points, double shade, int id){ std::vector<double> shades(points.size()); for (int i = 0; i<(int)points.size(); i++)shades[i] = shade; drawCloud(points, shades, id); }
  void drawLines(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends);
  void drawCylinders(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &radii, int id);
  void drawEllipsoids(const std::vector<Eigen::Vector3d> &centres, const std::vector<Eigen::Matrix3d> &poses, const std::vector<Eigen::Vector3d> &radii, const Eigen::Vector3d &colour, int id);
  std::string fixedFrameId_;

  #if defined(USE_ROS)
  ros::NodeHandle n;
  ros::Publisher cloudPublisher[2];
  ros::Publisher linePublisher;
  ros::Publisher cylinderPublisher[2];
  ros::Publisher ellipsoidPublisher[6];
  #endif
};
}
