// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYBRANCH_H
#define RAYLIB_RAYBRANCH_H

#include "../raycloud.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"

namespace ray
{
// Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
static const double minimum_score = 40.0;
static const double trunk_height_to_width = 4.0;  // height extent relative to real diameter of branch
static const double boundary_radius_scale = 3.0;  // how much farther out is the expected boundary compared to real
                                                  // branch radius? Larger requires more space to declare it a branch

struct Trunk
{
  Trunk();
  Eigen::Vector3d centre;
  double radius;
  double score;
  double last_score;
  double length;
  double actual_length;
  Eigen::Vector3d dir;
  int parent;
  double tree_score;
  double distance_to_ground;
  double ground_height;
  bool active;
  bool visited;

  void getOverlap(const Grid<Eigen::Vector3d> &grid, std::vector<Eigen::Vector3d> &points, double spacing);
  void estimatePose(const std::vector<Eigen::Vector3d> &points);
  void updateDirection(const std::vector<Eigen::Vector3d> &points);
  void updateCentre(const std::vector<Eigen::Vector3d> &points);
  void updateRadiusAndScore(const std::vector<Eigen::Vector3d> &points);
};

// grid structure as input to topology optimisation
struct RayGrid2D
{
  void init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width)
  {
    min_bound_ = min_bound;
    pixel_width_ = pixel_width;
    Eigen::Vector3d extent = max_bound - min_bound;
    dims_ = Eigen::Vector3d(std::ceil(extent[0] / pixel_width), std::ceil(extent[1] / pixel_width),
                            std::ceil(extent[2] / pixel_width))
              .cast<int>();
    std::cout << "min: " << min_bound.transpose() << ", ext: " << extent.transpose() << ", dims: " << dims_.transpose()
              << std::endl;

    pixels_.resize(dims_[0] * dims_[1]);
    memset(&pixels_[0], 0, sizeof(Pixel) * pixels_.size());
  }
  struct Pixel
  {
    std::vector<int> ray_ids;
    bool filled;
  };
  Eigen::Vector3i dims_;
  Eigen::Vector3d min_bound_;
  double pixel_width_;
  std::vector<Pixel> pixels_;
  Pixel dummy_pixel_;
  inline Eigen::Vector3i pixelIndex(const Eigen::Vector3d &pos) const
  {
    return ((pos - min_bound_) / pixel_width_).cast<int>();
  }
  inline const Pixel &pixel(const Eigen::Vector3i &index) const
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline Pixel &pixel(const Eigen::Vector3i &index)
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline const Pixel &pixel(const Eigen::Vector3d &pos) const { return pixel(pixelIndex(pos)); }
  inline Pixel &pixel(const Eigen::Vector3d &pos) { return pixel(pixelIndex(pos)); }

  // takes the filled cells and adds the ray ids that overlap these filled cells
  void fillRays(const Cloud &cloud)
  {
    ray::Cuboid bounds_;
    const double eps = 1e-9;
    bounds_.min_bound_ = min_bound_ + Eigen::Vector3d(eps, eps, eps);
    bounds_.max_bound_ = min_bound_ + dims_.cast<double>() * pixel_width_ - Eigen::Vector3d(eps, eps, eps);

    for (size_t i = 0; i < cloud.ends.size(); ++i)
    {
      Eigen::Vector3d start = cloud.starts[i];
      Eigen::Vector3d end = cloud.ends[i];
      bounds_.clipRay(start, end);

      // now walk the pixels
      const Eigen::Vector3d dir = end - start;
      const Eigen::Vector3d source = (start - min_bound_) / pixel_width_;
      const Eigen::Vector3d target = (end - min_bound_) / pixel_width_;
      const double length = dir.norm();
      const double eps = 1e-9;  // to stay away from edge cases
      const double maxDist =
        (target - source).norm() - 2.0;  // remove 2 subpixels to give a small buffer around the object

      // cached values to speed up the loop below
      Eigen::Vector3i adds;
      Eigen::Vector3d offsets;
      for (int k = 0; k < 3; ++k)
      {
        if (dir[k] > 0.0)
        {
          adds[k] = 1;
          offsets[k] = 0.5;
        }
        else
        {
          adds[k] = -1;
          offsets[k] = -0.5;
        }
      }

      Eigen::Vector3d p = source;  // our moving variable as we walk over the grid
      Eigen::Vector3i inds = p.cast<int>();
      double depth = 0;
      // walk over the grid, one pixel at a time.
      do
      {
        double ls[2] = { (round(p[0] + offsets[0]) - p[0]) / dir[0], (round(p[1] + offsets[1]) - p[1]) / dir[1] };
        int axis = (ls[0] < ls[1]) ? 0 : 1;
        inds[axis] += adds[axis];
        if (inds[axis] < 0 || inds[axis] >= dims_[axis])
        {
          break;
        }
        double minL = ls[axis] * length;
        depth += minL + eps;
        p = source + dir * (depth / length);
        Pixel &pix = pixel(inds);
        if (pix.filled)
        {
          pix.ray_ids.push_back((int)i);
        }
      } while (depth <= maxDist);
    }
  }
};

}  // namespace ray
#endif  // RAYLIB_RAYBRANCH_H