// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYGRID2D_H
#define RAYLIB_RAYGRID2D_H

#include "../rayutils.h"
#include "raylib/raylibconfig.h"
#include "raytrunks.h"
#define GRID2D_SUBPIXELS 4
namespace ray
{
// 2D occupancy grid structure that stores occupancy (density of pixel overlapping rays)
// If pixel sparsity is the number of subpixels that overlap rays, divided by the number of subpixels, then
// pixel density is 1 minus sparsity. A pixel with lots of rays passing through it is a strong indication that
// it is free space, so its density will be low.   
class RAYLIB_EXPORT OccupancyGrid2D
{
public:
  /// initialise for a given bounds and pixel width
  void init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width);

  /// save the grid
  void save(const std::string &filename);

  /// load the grid
  bool load(const std::string &filename);

  /// pixel structure that is used to estimate the density of free space
  struct Pixel
  {
    uint16_t bits;
    inline double density() const
    {
      return 1.0 - (static_cast<double>(bits) / static_cast<double>(GRID2D_SUBPIXELS * GRID2D_SUBPIXELS));
    }
  };

  /// 3D index of the pixel
  inline Eigen::Vector3i pixelIndex(const Eigen::Vector3d &pos) const
  {
    return ((pos - min_bound_) / pixel_width_).cast<int>();
  }
  /// pixel accessors
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

  /// fill in the occupancy data (pixels_) based on the rays in the cloud @c cloudname,
  /// within a height window @c clip_min to @c clip_max
  void fillDensities(const std::string &cloudname, const Eigen::ArrayXXd &lows, double clip_min, double clip_max);

  /// draw the occupancy grid
  void draw(const std::string &filename);

  const Eigen::Vector3i &dims() const { return dims_; }

private:
  Eigen::Vector3i dims_;
  Eigen::Vector3d min_bound_;
  double pixel_width_;
  std::vector<Pixel> pixels_;
  Pixel dummy_pixel_;
};

// A similar 2d grid structure, but this stores the ray indices per pixel
// this is an acceleration structure that allows overlapping rays to be quickly accessed at any location
class RAYLIB_EXPORT RayIndexGrid2D
{
public:
  void init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width);
  struct Pixel
  {
    std::vector<int> ray_ids;
    bool filled;
  };
  inline Eigen::Vector3i pixelIndex(const Eigen::Vector3d &pos) const
  {
    return ((pos - min_bound_) / pixel_width_).cast<int>();
  }
  inline const Pixel &pixel(const Eigen::Vector3i &index) const
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
    {
      return dummy_pixel_;
    }
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline Pixel &pixel(const Eigen::Vector3i &index)
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
    {
      return dummy_pixel_;
    }
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline const Pixel &pixel(const Eigen::Vector3d &pos) const { return pixel(pixelIndex(pos)); }
  inline Pixel &pixel(const Eigen::Vector3d &pos) { return pixel(pixelIndex(pos)); }

  // takes the filled cells and adds the ray ids that overlap these filled cells
  void fillRays(const Cloud &cloud);

private:
  Eigen::Vector3i dims_;       // dimensions of the grid. Only the first two elements are used here
  Eigen::Vector3d min_bound_;  // minimum bound of grid, SI units
  double pixel_width_;         // pixel width
  std::vector<Pixel> pixels_;  // storing the 2D data contiguously as a vector
  Pixel dummy_pixel_;          // to allow return values (containing no data) for out-of-range locations
};


}  // namespace ray

#endif  // RAYLIB_RAYGRID2D_H
