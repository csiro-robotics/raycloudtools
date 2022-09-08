// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYALIGNMENT_H
#define RAYLIB_RAYALIGNMENT_H

#include "raylib/raylibconfig.h"

#include "raycloud.h"
#include "rayutils.h"

#include <complex>

typedef std::complex<double> Complex;

namespace ray
{
/// Coarse raycloud alignment. This translates and 'yaw's the ray cloud, under the common assumption that pitch and roll
/// are already accurate. Transforms the first cloud in the pair @c clouds, to align with the second.
/// This is a cross-correlation method that requires a @c voxel_width (typically on the order of a metre)
/// the @c verbose argument saves out plan-view images at each step of the method.
/// The method uses a scale-free Fourier-Mellin transform to efficiently cross-correlate the cloud's end point
/// densities. NOTE @c clouds is a pair of clouds, it should point to an array with at least 2 elements
void RAYLIB_EXPORT alignCloud0ToCloud1(Cloud *clouds, double voxel_width, bool verbose = false);

/// 3D grid structure of complex numbers, for performing fast Fourier transforms (FFTs)
struct Array3D
{
  /// Initialise the grid with bounds, cell width and either a maximum bound or a dimensions vector
  void init(const Eigen::Vector3d &box_min, double voxel_width, const Eigen::Vector3i &dimensions);
  void init(const Eigen::Vector3d &box_min, const Eigen::Vector3d &box_max, double voxel_width);

  // Fast Fourier Transform
  void fft();
  // Inverse Fast Fourier Transform
  void inverseFft();

  void operator*=(const Array3D &other);

  // Accessors and modifiers
  inline Complex &operator()(int x, int y, int z) { return cells_[x + dims_[0] * y + dims_[0] * dims_[1] * z]; }
  inline const Complex &operator()(int x, int y, int z) const
  {
    return cells_[x + dims_[0] * y + dims_[0] * dims_[1] * z];
  }
  inline Complex &operator()(const Eigen::Vector3i &index) { return (*this)(index[0], index[1], index[2]); }
  inline const Complex &operator()(const Eigen::Vector3i &index) const { return (*this)(index[0], index[1], index[2]); }
  Complex &operator()(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d index = (pos - box_min_) / voxel_width_;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && index[0] < (double)dims_[0] &&
        index[1] < (double)dims_[1] && index[2] < (double)dims_[2])
      return (*this)(Eigen::Vector3i(index.cast<int>()));
    return null_cell_;
  }
  const Complex &operator()(const Eigen::Vector3d &pos) const
  {
    Eigen::Vector3d index = (pos - box_min_) / voxel_width_;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && index[0] < (double)dims_[0] &&
        index[1] < (double)dims_[1] && index[2] < (double)dims_[2])
      return (*this)(Eigen::Vector3i(index.cast<int>()));
    return null_cell_;
  }
  inline Eigen::Vector3i &dimensions() { return dims_; }
  inline const Eigen::Vector3i &dimensions() const { return dims_; }
  inline double voxelWidth() { return voxel_width_; }
  void clearCells() { cells_.clear(); }

  void conjugate();

  // Location in the grid of the cell wiith the largest real part
  Eigen::Vector3i maxRealIndex() const;

  // Fill grid based on the rays in the ray cloud
  void fillWithRays(const Cloud &cloud);

  Eigen::Vector3d box_min_, box_max_;
  double voxel_width_;

private:
  Eigen::Vector3i dims_;
  std::vector<Complex> cells_;
  Complex null_cell_;
};

}  // namespace ray

#endif  // RAYLIB_RAYALIGNMENT_H
