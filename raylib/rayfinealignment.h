// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYFINELIGNMENT_H
#define RAYLIB_RAYFINELIGNMENT_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"


namespace ray
{
#define state_size 12

class RAYLIB_EXPORT FineAlignment
{
public:
  FineAlignment(Cloud *clouds, bool non_rigid, bool verbose) : clouds_(clouds), non_rigid_(non_rigid), verbose_(verbose) {}
  void align();

private:
  struct RAYLIB_EXPORT Surfel
  {
    Surfel(){}
    Surfel(const Eigen::Vector3d &centroid, const Eigen::Matrix3d &matrix, const Eigen::Vector3d &width,
      const Eigen::Vector3d &normal, bool is_plane) :
        centroid(centroid), matrix(matrix), width(width), normal(normal), is_plane(is_plane) {}

    Eigen::Vector3d centroid;
    Eigen::Matrix3d matrix;
    Eigen::Vector3d width;
    Eigen::Vector3d normal;
    bool is_plane;
    static void draw(const std::vector<Surfel> &surfels, const Eigen::Vector3d &colour);
  };

  struct Match
  {
    int ids[2];
    Eigen::Vector3d normal;
  };

  struct QuadraticTransformation
  {
    QuadraticTransformation(){}
    QuadraticTransformation(const Eigen::Matrix<double, state_size, 1> &x)
    {
      translation = Eigen::Vector3d(x[0], x[1], x[2]);
      rotation = Eigen::Vector3d(x[3], x[4], x[5]);
      a = Eigen::Vector3d(x[6], x[7], 0);
      b = Eigen::Vector3d(x[8], x[9], 0);
      c = Eigen::Vector3d(x[10], x[11], 0);    
    }  
    Eigen::Vector3d translation, rotation; // the linear components
    Eigen::Vector3d a, b, c; // the principle quadratic components
    Pose getEuclideanPart() const 
    { 
      return Pose(translation, Eigen::Quaterniond(Eigen::AngleAxisd(rotation.norm(), rotation.normalized()))); 
    }
  };

  struct LinearSystem
  {
    LinearSystem(){ At_A.setZero(); At_b.setZero(); }
    Eigen::Matrix<double, state_size, 1> solve();
    Eigen::Matrix<double, state_size, state_size> At_A;
    Eigen::Matrix<double, state_size, 1> At_b;
  };

  void generateSurfels();
  void generateSurfelMatches(std::vector<Match> &matches);
  void buildLinearSystem(const std::vector<Match> &matches, double d, FineAlignment::LinearSystem &system);
  void updateLinearSystem(std::vector<Match> &matches, const QuadraticTransformation &trans);

  // primary data:
  Cloud *clouds_;
  double non_rigid_;
  double verbose_;
  const double max_normal_difference_ = 0.5;

  // derived data
  std::vector<Surfel> surfels_[2];
  double translation_weight_;
  Eigen::Vector3d centres_[2];
};
}  // namespace ray

#endif  // RAYLIB_RAYFINELIGNMENT_H
