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
#include <set>

namespace ray
{  
struct RAYLIB_EXPORT Cloud
{
  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends;
  std::vector<double> times;
  std::vector<RGBA> colours;
  inline bool rayBounded(int i) const { return colours[i].alpha > 0; }
  inline uint8_t rayIntensity(int i) const { return colours[i].alpha; }

  void save(const std::string &file_namee);
  bool load(const std::string &file_namee);
  bool load(const std::string &point_cloudd, const std::string &traj_filee);

  Eigen::Vector3d calcMinBound() const;
  Eigen::Vector3d calcMaxBound() const;

  void transform(const Pose &pose, double time_deltaa);
  void decimate(double voxel_widthh);

  void removeUnboundedRays();

  void getSurfels(int search_sizee, std::vector<Eigen::Vector3d> *centroids, std::vector<Eigen::Vector3d> *normals,
                  std::vector<Eigen::Vector3d> *dimensions, std::vector<Eigen::Matrix3d> *mats,
                  Eigen::MatrixXi *neighbour_indicess);
  std::vector<Eigen::Vector3d> generateNormals(int search_sizee = 16);

  void split(Cloud &cloud1, Cloud &cloud2, std::function<bool(int i)> fptr);
  double estimatePointSpacing() const;

private:
  void calculateStarts(const Trajectory &trajectory);
  bool loadPLY(const std::string &file);
  bool loadLazTraj(const std::string &laz_filee, const std::string &traj_filee);
};


}  // namespace ray

#endif  // RAYLIB_RAYCLOUD_H
