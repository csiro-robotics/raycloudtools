// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include "rayutils.h"
#include "raycloud.h"
#include <complex>
typedef std::complex<double> Complex;

namespace RAY
{
struct Array3D
{
  void init(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth);

  void FFT();
  void inverseFFT();

  void operator *=(const Array3D &other);

  inline Complex &operator()(const int &x, const int &y, const int &z)
  {
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }
  inline const Complex &operator()(const int &x, const int &y, const int &z) const
  {
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }
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
  inline Complex &operator()(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d index = (pos - boxMin)/voxelWidth;
    return (*this)(index[0], index[1], index[2]);
  }  
  inline const Complex &operator()(const Eigen::Vector3d &pos) const
  {
    Eigen::Vector3d index = (pos - boxMin)/voxelWidth;
    return (*this)(index[0], index[1], index[2]);
  }
  void conjugate();
  Eigen::Vector3i maxRealIndex() const;

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
  Eigen::Vector3i dims;
  std::vector<Complex> cells;
};

struct Array1D
{
  void init(int length);

  void FFT();
  void inverseFFT();

  void operator *=(const Array1D &other);
  inline Complex &operator()(const int &x){ return cells[x]; }
  inline const Complex &operator()(const int &x) const { return cells[x]; }
  inline void operator +=(const Array1D &other)
  {
    for (int i = 0; i<(int)cells.size(); i++)
      cells[i] += other.cells[i];
  }

  int maxRealIndex() const;
  void conjugate();

  std::vector<Complex> cells;
};

struct AlignTranslationYaw
{
  void alignCloud0ToCloud1(double voxelWidth);
  Cloud clouds[2];
};

}
