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
class Array3D
{
public:
  void init(const Eigen::Vector3d &boxMin, const Eigen::Vector3d &boxMax, double voxelWidth);

  Array3D &cwiseProductInplace(const Array3D &other);
  Array3D &conjugateInplace();
  Array3D &fft();
  Array3D &ifft();
  Array3D &convolve(Array3D &target);
  Eigen::Vector3i maxRealIndex() const;

  inline Complex &operator()(const int &x, const int &y, const int &z)
  {
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }
  /// Access the data for the requested co-ordinate.
  inline const Complex &operator()(const int &x, const int &y, const int &z) const
  {
    return cells[x + dims[0]*y + dims[0]*dims[1]*z];
  }

  inline Complex &operator()(const Eigen::Vector3i &index)
  {
//    ASSERT(x<0 || x>=dims[0] || y<0 || y>= dims[1] || z<0 || z>=dims[2]);
    return (*this)(index[0], index[1], index[2]);
  }
  /// Access the data for the requested co-ordinate.
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

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
  Eigen::Vector3i dims;
  std::vector<Complex> cells;
};

class Array1D
{
public:
  void init(int length);

  Array1D &cwiseProductInplace(const Array1D &other);
  Array1D &conjugateInplace();
  Array1D &fft();
  Array1D &ifft();
  int maxRealIndex() const;

  inline Complex &operator()(const int &x)
  {
    return cells[x];
  }
  inline const Complex &operator()(const int &x) const
  {
    return cells[x];
  }
  inline void operator +=(const Array1D &other)
  {
    for (int i = 0; i<(int)cells.size(); i++)
      cells[i] += other.cells[i];
  }

  std::vector<Complex> cells;
};



}
