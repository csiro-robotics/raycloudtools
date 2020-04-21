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
  Complex &operator()(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d index = (pos - boxMin)/voxelWidth;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && 
        index[0] < (double)dims[0] && index[1] < (double)dims[1] && index[2] < (double)dims[2])
      return (*this)(index[0], index[1], index[2]);
    return nullCell;
  }  
  const Complex &operator()(const Eigen::Vector3d &pos) const
  {
    Eigen::Vector3d index = (pos - boxMin)/voxelWidth;
    if (index[0] >= 0.0 && index[1] >= 0.0 && index[2] >= 0.0 && 
        index[0] < (double)dims[0] && index[1] < (double)dims[1] && index[2] < (double)dims[2])
      return (*this)(index[0], index[1], index[2]);
    return nullCell;
  }
  inline void stableDivideBy(const Array3D &other, double factor)
  {
 /*   double avg = 0.0, avg2 = 0;
    double count = 0.0;
    int numNeg = 0;
    double maxVal = 0.0;
    for (int i = 0; i<(int)cells.size(); i++)
    {
      double w = other.cells[i].real(); 
      if (w < 0.0)
        numNeg++;
      maxVal = std::max(maxVal, w);
      if (w > 0)
      {
        avg += w;
        avg2 += w*w;
        count++;
      }
    }
    avg2 /= avg;
    avg /= count;
    std::cout << "factor supplied: " << factor << ", max weight: " << maxVal << ", average weight: " << avg << ", average weigthed weight: " << avg2 << ", num negative weights: " << numNeg << "/" << cells.size() << std::endl;
  */
    for (int i = 0; i<(int)cells.size(); i++)
      cells[i] /= factor + std::max(0.0, other.cells[i].real()); 
  }
  void conjugate();
  Eigen::Vector3i maxRealIndex() const;
  void fillWithRays(const Cloud &cloud);

  Eigen::Vector3d boxMin, boxMax;
  double voxelWidth;
  Eigen::Vector3i dims;
  std::vector<Complex> cells;
  Complex nullCell;
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
  inline void stableDivideBy(const Array1D &other, double factor)
  {
    for (int i = 0; i<(int)cells.size(); i++)
      cells[i] /= factor + std::max(0.0, other.cells[i].real()); 
  }
  void polarCrossCorrelation(const Array3D *arrays, const Array3D *weights, bool verbose);

  int maxRealIndex() const;
  void conjugate();

  std::vector<Complex> cells;
};

struct AlignTranslationYaw
{
  void alignCloud0ToCloud1(double voxelWidth, bool verbose=false);
  Cloud clouds[2];
};

}
