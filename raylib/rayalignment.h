// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYALIGNMENT_H
#define RAYLIB_RAYALIGNMENT_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raycloud.h"

#include <complex>

typedef std::complex<double> Complex;

namespace ray
{
class RAYLIB_EXPORT Array3D
{
public:
  void init(const Eigen::Vector3d &box_minn, const Eigen::Vector3d &box_maxx, double voxel_widthh);

  void fft();
  void inverseFft();

  void operator*=(const Array3D &other);

  inline Complex &operator()(int x, int y, int z) { return cells[x + dims[0] * y + dims[0] * dims[1] * z]; }
  inline const Complex &operator()(int x, int y, int z) const { return cells[x + dims[0] * y + dims[0] * dims[1] * z]; }
  inline Complex &operator()(const Eigen::Vector3i &index)
  {
    //    ASSERT(x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2]);
    return (*this)(index[0], index[1], index[2]);
  }
  inline const Complex &operator()(const Eigen::Vector3i &index) const
  {
    //    ASSERT(x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2]);
    return (*this)(index[0], index[1], index[2]);
  }
  Complex &operator()(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d index = (pos - box_min) / voxel_width;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && index[0] < (double)dims[0] &&
        index[1] < (double)dims[1] && index[2] < (double)dims[2])
      return (*this)(Eigen::Vector3i(index.cast<int>()));
    return null_cell;
  }
  const Complex &operator()(const Eigen::Vector3d &pos) const
  {
    Eigen::Vector3d index = (pos - box_min) / voxel_width;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && index[0] < (double)dims[0] &&
        index[1] < (double)dims[1] && index[2] < (double)dims[2])
      return (*this)(Eigen::Vector3i(index.cast<int>()));
    return null_cell;
  }
  void conjugate();
  Eigen::Vector3i maxRealIndex() const;
  void fillWithRays(const Cloud &cloud);

  Eigen::Vector3d box_min, box_max;
  double voxel_width;
  Eigen::Vector3i dims;
  std::vector<Complex> cells;
  Complex null_cell;
};

class RAYLIB_EXPORT Array1D
{
public:
  void init(int length);

  void fft();
  void inverseFft();

  void operator*=(const Array1D &other);
  inline Complex &operator()(const int &x) { return cells[x]; }
  inline const Complex &operator()(const int &x) const { return cells[x]; }
  inline void operator+=(const Array1D &other)
  {
    for (int i = 0; i < (int)cells.size(); i++) cells[i] += other.cells[i];
  }
  void polarCrossCorrelation(const Array3D *arrays, bool verbose);

  int maxRealIndex() const;
  void conjugate();

  std::vector<Complex> cells;
};

class AlignTranslationYaw
{
public:
  void alignCloud0ToCloud1(double voxel_widthh, bool verbose = false);
  Cloud clouds[2];
};

}  // namespace ray

#endif  // RAYLIB_RAYALIGNMENT_H
