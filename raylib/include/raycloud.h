// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"
#include "raygrid.h"

namespace RAY
{
struct Ellipsoid
{
  Eigen::Vector3d pos;
  Eigen::Vector3d vectors[3];
  double time;
  double size;
  bool transient;
};

struct Cloud
{
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends; 
  std::vector<double> times;
  std::vector<RGBA> colours; 
  inline bool rayBounded(int i){ return colours[i].alpha > 0; }
  inline uint8_t rayIntensity(int i){ return colours[i].alpha; }

  void save(const std::string &fileName);
  bool load(const std::string &fileName);
  bool load(const std::string &pointCloud, const std::string &trajFile);

  void transform(const Pose &pose, double timeDelta);
  void decimate(double voxelWidth);

  void removeUnboundedRays();
  std::vector<Eigen::Vector3d> generateNormals(int searchSize = 16);
  void findTransients(Cloud &transient, Cloud &fixed, double timeDelta, const std::string &mergeType);
  void combine(std::vector<Cloud> &clouds, Cloud &differences, const std::string &mergeType);
  void markIntersectedEllipsoids(Grid<Ellipsoid *> &grid, double timeDelta, const std::string &mergeType);
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids);

  Eigen::Vector3d calcMinBound();
  Eigen::Vector3d calcMaxBound();

protected:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &lazFile, const std::string &trajFile);
};


}
