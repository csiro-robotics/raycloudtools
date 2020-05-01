// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUD_H
#define RAYLIB_RAYCLOUD_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raypose.h"
#include "raytrajectory.h"
#include "raygrid.h"

namespace ray
{
struct RAYLIB_EXPORT Ellipsoid
{
  Eigen::Vector3d pos;
  Eigen::Vector3d vectors[3];
  double time;
  Eigen::Vector3d extents;
  double opacity;
  double planarity;
  size_t num_rays;
  size_t num_gone;
  inline void setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals)
  {
    double max_rr = std::max(vals[0], std::max(vals[1], vals[2]));
    const Eigen::Vector3d &x = vecs.col(0);
    const Eigen::Vector3d &y = vecs.col(1);
    const Eigen::Vector3d &z = vecs.col(2);
    extents[0] = std::min(max_rr, abs(x[0])*vals[0] + abs(y[0])*vals[1] + abs(z[0])*vals[2]);
    extents[1] = std::min(max_rr, abs(x[1])*vals[0] + abs(y[1])*vals[1] + abs(z[1])*vals[2]);
    extents[2] = std::min(max_rr, abs(x[2])*vals[0] + abs(y[2])*vals[1] + abs(z[2])*vals[2]);
  }
  void setPlanarity(const Eigen::Vector3d &vals)
  {
    planarity = (vals[1]-vals[0])/vals[1];
  }
  bool transient;
};

struct RAYLIB_EXPORT Cloud
{
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends; 
  std::vector<double> times;
  std::vector<RGBA> colours; 
  inline bool rayBounded(int i){ return colours[i].alpha > 0; }
  inline uint8_t rayIntensity(int i){ return colours[i].alpha; }

  void save(const std::string &file_namee);
  bool load(const std::string &file_namee);
  bool load(const std::string &point_cloudd, const std::string &traj_filee);

  void transform(const Pose &pose, double time_deltaa);
  void decimate(double voxel_widthh);

  void removeUnboundedRays();
  std::vector<Eigen::Vector3d> generateNormals(int search_sizee = 16);
  void findTransients(Cloud &transient, Cloud &fixed, const std::string &merge_typee, double num_rays, bool colour_cloudd);
  void combine(std::vector<Cloud> &clouds, Cloud &differences, const std::string &merge_typee, double num_rays);
  void markIntersectedEllipsoids(Grid<int> &grid, std::vector<bool> &transients, std::vector<Ellipsoid> &ellipsoids, const std::string &merge_typee, double num_rays, bool self_transientt);
  void generateEllipsoids(std::vector<Ellipsoid> &ellipsoids);
  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);
  void getSurfels(int search_sizee, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals, std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats, Eigen::MatrixXi *neighbour_indicess);

  Eigen::Vector3d calcMinBound();
  Eigen::Vector3d calcMaxBound();

private:  
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &laz_filee, const std::string &traj_filee);
};


}

#endif // RAYLIB_RAYCLOUD_H
