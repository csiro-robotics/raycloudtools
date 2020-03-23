// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayalignment.h"
#include "simple_fft/fft.h"

using namespace std;
using namespace Eigen;
using namespace RAY;

void Array3D::init(const Vector3d &boxMin, const Vector3d &boxMax, double voxelWidth)
{
  this->boxMin = boxMin;
  this->boxMax = boxMax;
  this->voxelWidth = voxelWidth;
  Eigen::Vector3d diff = (boxMax - boxMin)/voxelWidth;
  // HERE we need to make it a power of two
  for (int i = 0; i<3; i++)
    dims[i] = 1 << (int)ceil(log2(ceil(diff[i]))); // next power of two larger than diff
  cells.resize(dims[0]*dims[1]*dims[2]);
  memset(&cells[0], 0, cells.size()*sizeof(Complex));
}

Array3D &Array3D::cwiseProductInplace(const Array3D &other)
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] *= other.cells[i];
  return *this;
}

Array3D &Array3D::conjugateInplace()
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] = conj(cells[i]);
  return *this;
}

Array3D &Array3D::fft()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate FFT: " << error << endl;
  return *this;
}

Array3D &Array3D::ifft()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate inverse FFT: " << error << endl;
  return *this;
}

Array3D &Array3D::convolve(Array3D &target)
{
  return fft().cwiseProductInplace(target.fft().conjugateInplace()).ifft();
}

Vector3i Array3D::maxRealIndex() const
{
  Vector3i index;
  double highest = numeric_limits<double>::lowest();
  for (int i = 0; i<(int)cells.size(); i++)
  {
    const double &score = cells[i].real();
    if (score > highest)
    {
      index = Vector3i(i%dims[0], (i/dims[0])%dims[1], i/(dims[0]*dims[1]));
      highest = score;
    }
  }
  return index;
}

/**************************************************************************************************/

void Array1D::init(int length)
{
  cells.resize(length);
  memset(&cells[0], 0, cells.size()*sizeof(Complex));
}

Array1D &Array1D::cwiseProductInplace(const Array1D &other)
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] *= other.cells[i];
  return *this;
}

Array1D &Array1D::conjugateInplace()
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] = conj(cells[i]);
  return *this;
}

Array1D &Array1D::fft()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, cells.size(), error))
    cout << "failed to calculate FFT: " << error << endl;
  return *this;
}

Array1D &Array1D::ifft()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, cells.size(), error))
    cout << "failed to calculate inverse FFT: " << error << endl;
  return *this;
}


int Array1D::maxRealIndex() const
{
  int index;
  double highest = numeric_limits<double>::lowest();
  for (int i = 0; i<(int)cells.size(); i++)
  {
    const double &score = cells[i].real();
    if (score > highest)
    {
      index = i;
      highest = score;
    }
  }
  return index;
}
