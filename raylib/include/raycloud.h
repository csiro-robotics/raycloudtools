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
  Eigen::Vector3d extents;
  double opacity;
  double planarity;
  int numRays;
  int numGone;
  inline void setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals)
  {
    double maxR = std::max(vals[0], std::max(vals[1], vals[2]));
    const Eigen::Vector3d &x = vecs.col(0);
    const Eigen::Vector3d &y = vecs.col(1);
    const Eigen::Vector3d &z = vecs.col(2);
    extents[0] = std::min(maxR, abs(x[0])*vals[0] + abs(y[0])*vals[1] + abs(z[0])*vals[2]);
    extents[1] = std::min(maxR, abs(x[1])*vals[0] + abs(y[1])*vals[1] + abs(z[1])*vals[2]);
    extents[2] = std::min(maxR, abs(x[2])*vals[0] + abs(y[2])*vals[1] + abs(z[2])*vals[2]);
  }
  void setPlanarity(const Eigen::Vector3d &vals)
  {
    planarity = (vals[1]-vals[0])/vals[1];
  }
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
  void findTransients(Cloud &transient, Cloud &fixed, const std::string &mergeType, double numRays, bool colourCloud);
  void combine(std::vector<Cloud> &clouds, Cloud &differences, const std::string &mergeType, double numRays);
  void markIntersectedEllipsoids(Grid<int> &grid, std::vector<bool> &transients, std::vector<Ellipsoid> &ellipsoids, const std::string &mergeType, double numRays, bool selfTransient);
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids);
  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);
  void getSurfels(int searchSize, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals, std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats, Eigen::MatrixXi *neighbourIndices);

  Eigen::Vector3d calcMinBound();
  Eigen::Vector3d calcMaxBound();

private:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &lazFile, const std::string &trajFile);
};


}
