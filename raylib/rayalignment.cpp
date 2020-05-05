// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayalignment.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "imagewrite.h"

#include "simple_fft/fft.h"

#include <cinttypes>
#include <iostream>
#include <complex>

using namespace std;
using namespace Eigen;
using namespace ray;

typedef complex<double> Complex;
static const double kHighPassPower = 0.25;  // This fixes inout->inout11, inoutD->inoutB2 and house_inside->house3.
                                            // Doesn't break any. power=0.25. 0 is turned off.

struct Col
{
  uint8_t r, g, b, a;
};

void Array3D::init(const Vector3d &box_min, const Vector3d &box_max, double voxel_width)
{
  this->box_min = box_min;
  this->box_max = box_max;
  this->voxel_width = voxel_width;
  Eigen::Vector3d diff = (box_max - box_min) / voxel_width;
  // HERE we need to make it a power of two
  for (int i = 0; i < 3; i++) dims[i] = 1 << (int)ceil(log2(ceil(diff[i])));  // next power of two larger than diff
  cells.resize(dims[0] * dims[1] * dims[2]);
  memset(&cells[0], 0, cells.size() * sizeof(Complex));
  null_cell = 0;
}

void Array3D::operator*=(const Array3D &other)
{
  for (int i = 0; i < (int)cells.size(); i++) cells[i] *= other.cells[i];
}

void Array3D::conjugate()
{
  for (int i = 0; i < (int)cells.size(); i++) cells[i] = conj(cells[i]);
}

void Array3D::fft()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate FFT: " << error << endl;
}

void Array3D::inverseFft()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate inverse FFT: " << error << endl;
}

Vector3i Array3D::maxRealIndex() const
{
  Vector3i index;
  double highest = numeric_limits<double>::lowest();
  for (int i = 0; i < (int)cells.size(); i++)
  {
    const double &score = cells[i].real();
    if (score > highest)
    {
      index = Vector3i(i % dims[0], (i / dims[0]) % dims[1], i / (dims[0] * dims[1]));
      highest = score;
    }
  }
  return index;
}

void Array3D::fillWithRays(const Cloud &cloud)
{
  // unlike the end point densities, the weight is just 0 or 1, but requires walking through the grid for every ray
  // maybe a better choice would be a reuseable 'volume' function (occupancy grid).
  for (int i = 0; i < (int)cloud.ends.size(); i++)
  {
    // TODO: Should be a lambda function!
    Vector3d dir = cloud.ends[i] - cloud.starts[i];
    Vector3d dir_sign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (cloud.starts[i] - box_min) / voxel_width;
    Vector3d end = (cloud.ends[i] - box_min) / voxel_width;
    Vector3i start_index(start.cast<int>());
    Vector3i end_index(end.cast<int>());
    double length_sqr = (end_index - start_index).squaredNorm();
    Vector3i index = start_index;
    while ((index - start_index).squaredNorm() <= length_sqr + 1e-10)
    {
      if (index[0] >= 0 && index[0] < dims[0] && index[1] >= 0 && index[1] < dims[1] && index[2] >= 0 &&
          index[2] < dims[2])
        (*this)(index[0], index[1], index[2]) += Complex(1, 0);  // add weight to these areas...

      Vector3d mid = box_min + voxel_width * Vector3d(index[0] + 0.5, index[1] + 0.5, index[2] + 0.5);
      Vector3d next_boundary = mid + 0.5 * voxel_width * dir_sign;
      Vector3d delta = next_boundary - cloud.starts[i];
      Vector3d d(delta[0] / dir[0], delta[1] / dir[1], delta[2] / dir[2]);
      if (d[0] < d[1] && d[0] < d[2])
        index[0] += dir_sign.cast<int>()[0];
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += dir_sign.cast<int>()[1];
      else
        index[2] += dir_sign.cast<int>()[2];
    }
  }
}

/**************************************************************************************************/

void Array1D::init(int length)
{
  cells.resize(length);
  memset(&cells[0], 0, cells.size() * sizeof(Complex));
}

void Array1D::operator*=(const Array1D &other)
{
  for (int i = 0; i < (int)cells.size(); i++) cells[i] *= other.cells[i];
}

void Array1D::conjugate()
{
  for (int i = 0; i < (int)cells.size(); i++) cells[i] = conj(cells[i]);
}

void Array1D::fft()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, cells.size(), error))
    cout << "failed to calculate FFT: " << error << endl;
}

void Array1D::inverseFft()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, cells.size(), error))
    cout << "failed to calculate inverse FFT: " << error << endl;
}

int Array1D::maxRealIndex() const
{
  int index = 0;
  double highest = numeric_limits<double>::lowest();
  for (int i = 0; i < (int)cells.size(); i++)
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

void drawArray(const Array3D &array, const Vector3i &dims, const string &file_name, int index)
{
  int width = dims[0];
  int height = dims[1];
  double max_val = 0.0;
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      double val = 0.0;
      for (int z = 0; z < dims[2]; z++) val += abs(array(x, y, z));
      max_val = max(max_val, val);
    }
  }

  vector<Col> pixels(width * height);
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Vector3d colour(0, 0, 0);
      for (int z = 0; z < dims[2]; z++)
      {
        double h = (double)z / (double)dims[2];
        Vector3d col;
        col[0] = 1.0 - h;
        col[2] = h;
        col[1] = 3.0 * col[0] * col[2];
        colour += abs(array(x, y, z)) * col;
      }
      colour *= 15.0 * 255.0 / max_val;
      Col col;
      col.r = uint8_t(clamped((int)colour[0], 0, 255));
      col.g = uint8_t(clamped((int)colour[1], 0, 255));
      col.b = uint8_t(clamped((int)colour[2], 0, 255));
      col.a = 255;
      pixels[(x + width / 2) % width + width * ((y + height / 2) % height)] = col;
    }
  }
  stringstream str;
  str << file_name << index << ".png";
  stbi_write_png(str.str().c_str(), width, height, 4, (void *)&pixels[0], 4 * width);
}

void drawArray(const vector<Array1D> &arrays, const Vector3i &dims, const string &file_name, int index)
{
  int width = dims[0];
  int height = dims[1];
  double max_val = 0.0;
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      double val = 0.0;
      for (int z = 0; z < dims[2]; z++) val += abs(arrays[y + dims[1] * z](x));
      max_val = max(max_val, val);
    }
  }

  vector<Col> pixels(width * height);
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Vector3d colour(0, 0, 0);
      for (int z = 0; z < dims[2]; z++)
      {
        double h = (double)z / (double)dims[2];
        Vector3d col;
        col[0] = 1.0 - h;
        col[2] = h;
        col[1] = 3.0 * col[0] * col[2];
        colour += abs(arrays[y + dims[1] * z](x)) * col;
      }
      colour *= 3.0 * 255.0 / max_val;
      Col col;
      col.r = uint8_t(clamped((int)colour[0], 0, 255));
      col.g = uint8_t(clamped((int)colour[1], 0, 255));
      col.b = uint8_t(clamped((int)colour[2], 0, 255));
      col.a = 255;
      pixels[(x + width / 2) % width + width * y] = col;
    }
  }
  stringstream str;
  str << file_name << index << ".png";
  stbi_write_png(str.str().c_str(), width, height, 4, (void *)&pixels[0], 4 * width);
}

void Array1D::polarCrossCorrelation(const Array3D *arrays, bool verbose)
{
  // OK cool, so next I need to re-map the two arrays into 4x1 grids...
  int max_rad = max(arrays[0].dims[0], arrays[0].dims[1]) / 2;
  Vector3i polar_dims = Vector3i(4 * max_rad, max_rad, arrays[0].dims[2]);
  vector<Array1D> polars[2];
  for (int c = 0; c < 2; c++)
  {
    vector<Array1D> &polar = polars[c];
    const Array3D &a = arrays[c];
    polar.resize(polar_dims[1] * polar_dims[2]);
    for (int j = 0; j < polar_dims[1]; j++)
      for (int k = 0; k < polar_dims[2]; k++) polar[j + polar_dims[1] * k].init(polar_dims[0]);

    // now map...
    for (int i = 0; i < polar_dims[0]; i++)
    {
      double angle = 2.0 * kPi * (double)(i + 0.5) / (double)polar_dims[0];
      for (int j = 0; j < polar_dims[1]; j++)
      {
        double radius = (0.5 + (double)j) / (double)polar_dims[1];
        Vector2d pos = radius * 0.5 * Vector2d((double)a.dims[0] * sin(angle), (double)a.dims[1] * cos(angle));
        if (pos[0] < 0.0)
          pos[0] += a.dims[0];
        if (pos[1] < 0.0)
          pos[1] += a.dims[1];
        int x = pos.cast<int>()[0];
        int y = pos.cast<int>()[1];
        int x2 = (x + 1) % a.dims[0];
        int y2 = (y + 1) % a.dims[1];
        double blend_x = pos[0] - (double)x;
        double blend_y = pos[1] - (double)y;
        for (int z = 0; z < polar_dims[2]; z++)
        {
          // bilinear interpolation -- for some reason LERP after abs is better than before abs
          double val = abs(a(x, y, z)) * (1.0 - blend_x) * (1.0 - blend_y) +
                       abs(a(x2, y, z)) * blend_x * (1.0 - blend_y) + abs(a(x, y2, z)) * (1.0 - blend_x) * blend_y +
                       abs(a(x2, y2, z)) * blend_x * blend_y;
          polar[j + polar_dims[1] * z](i) = Complex(radius * val, 0);
        }
      }
    }
    if (verbose)
      drawArray(polar, polar_dims, "translationInvPolar", c);
    for (int j = 0; j < polar_dims[1]; j++)
    {
      for (int k = 0; k < polar_dims[2]; k++)
      {
        int i = j + polar_dims[1] * k;
        polar[i].fft();
        if (kHighPassPower > 0.0)
        {
          for (int l = 0; l < (int)polar[i].cells.size(); l++)
            polar[i].cells[l] *= pow(min((double)l, (double)(polar[i].cells.size() - l)), kHighPassPower);
        }
      }
    }
    if (verbose)
      drawArray(polar, polar_dims, "euclideanInvariant", c);
  }

  // now get the inverse fft in place:
  init(polar_dims[0]);
  for (int i = 0; i < (int)polars[0].size(); i++)
  {
    polars[1][i].conjugate();
    polars[0][i] *= polars[1][i];
    polars[0][i].inverseFft();
    (*this) += polars[0][i];  // add all the results together into the first array
  }
}

/************************************************************************************/

void AlignTranslationYaw::alignCloud0ToCloud1(double voxel_width, bool verbose)
{
  // first we need to decimate the clouds into intensity grids..
  // I need to get a maximum box width, and individual box_min, boxMaxs
  Vector3d box_mins[2], box_width(0, 0, 0);
  for (int c = 0; c < 2; c++)
  {
    Vector3d box_min(1e10, 1e10, 1e10), box_max(-1e10, -1e10, -1e10);
    for (int i = 0; i < (int)clouds[c].ends.size(); i++)
    {
      if (clouds[c].rayBounded(i))
      {
        box_min = minVector(box_min, clouds[c].ends[i]);
        box_max = maxVector(box_max, clouds[c].ends[i]);
      }
    }
    box_mins[c] = box_min;
    Vector3d width = box_max - box_min;
    box_width = maxVector(box_width, width);
  }

  bool rotation_to_estimate = true;  // If we know there is no rotation between the clouds then we can save some cost

  Array3D arrays[2];
  // Now fill in the arrays with point density
  for (int c = 0; c < 2; c++)
  {
    arrays[c].init(box_mins[c], box_mins[c] + box_width, voxel_width);
    for (int i = 0; i < (int)clouds[c].ends.size(); i++)
      if (clouds[c].rayBounded(i))
        arrays[c](clouds[c].ends[i]) += Complex(1, 0);
    arrays[c].fft();
    if (verbose)
      drawArray(arrays[c], arrays[c].dims, "translationInvariant", c);
  }

  if (rotation_to_estimate)
  {
    Array1D polar;
    polar.polarCrossCorrelation(arrays, verbose);

    // get the angle of rotation
    int index = polar.maxRealIndex();
    // add a little bit of sub-pixel accuracy:
    double angle;
    int dim = int(polar.cells.size());
    int back = (index + dim - 1) % dim;
    int fwd = (index + 1) % dim;
    double y0 = polar(back).real();
    double y1 = polar(index).real();
    double y2 = polar(fwd).real();
    angle = index + 0.5 * (y0 - y2) / (y0 + y2 - 2.0 * y1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
    // but the FFT wraps around, so:
    if (angle > dim / 2)
      angle -= dim;
    angle *= 2.0 * kPi / (double)polar.cells.size();
    cout << "estimated yaw rotation: " << angle << endl;

    // ok, so let's rotate A towards B, and re-run the translation FFT
    clouds[0].transform(Pose(Vector3d(0, 0, 0), Quaterniond(AngleAxisd(angle, Vector3d(0, 0, 1)))), 0.0);

    box_mins[0] = Vector3d(1e10, 1e10, 1e10);
    for (int i = 0; i < (int)clouds[0].ends.size(); i++)
      if (clouds[0].rayBounded(i))
        box_mins[0] = minVector(box_mins[0], clouds[0].ends[i]);
    arrays[0].cells.clear();
    arrays[0].init(box_mins[0], box_mins[0] + box_width, voxel_width);

    for (int i = 0; i < (int)clouds[0].ends.size(); i++)
      if (clouds[0].rayBounded(i))
        arrays[0](clouds[0].ends[i]) += Complex(1, 0);

    arrays[0].fft();
    if (verbose)
      drawArray(arrays[0], arrays[0].dims, "translationInvariantWeighted", 0);
  }

  if (kHighPassPower > 0.0)
  {
    for (int c = 0; c < 2; c++)
    {
      for (int x = 0; x < arrays[c].dims[0]; x++)
      {
        double coord_x = x < arrays[c].dims[0] / 2 ? x : arrays[c].dims[0] - x;
        for (int y = 0; y < arrays[c].dims[1]; y++)
        {
          double coord_y = y < arrays[c].dims[1] / 2 ? y : arrays[c].dims[1] - y;
          for (int z = 0; z < arrays[c].dims[2]; z++)
          {
            double coord_z = z < arrays[c].dims[2] / 2 ? z : arrays[c].dims[2] - z;
            arrays[c](x, y, z) *= pow(sqr(coord_x) + sqr(coord_y) + sqr(coord_z), kHighPassPower);
          }
        }
      }
      if (verbose)
        drawArray(arrays[c], arrays[c].dims, "normalised", c);
    }
  }
  /****************************************************************************************************/
  // now get the the translation part
  arrays[1].conjugate();
  arrays[0] *= arrays[1];
  arrays[0].inverseFft();

  // find the peak
  Array3D &array = arrays[0];
  Vector3i ind = array.maxRealIndex();
  // add a little bit of sub-pixel accuracy:
  Vector3d pos;
  for (int axis = 0; axis < 3; axis++)
  {
    Vector3i back = ind, fwd = ind;
    int &dim = array.dims[axis];
    back[axis] = (ind[axis] + dim - 1) % dim;
    fwd[axis] = (ind[axis] + 1) % dim;
    double y0 = array(back).real();
    double y1 = array(ind).real();
    double y2 = array(fwd).real();
    pos[axis] =
      ind[axis] + 0.5 * (y0 - y2) / (y0 + y2 - 2.0 * y1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
    // but the FFT wraps around, so:
    if (pos[axis] > dim / 2)
      pos[axis] -= dim;
  }
  pos *= -array.voxel_width;
  pos += box_mins[1] - box_mins[0];
  cout << "estimated translation: " << pos.transpose() << endl;

  Pose transform(pos, Quaterniond::Identity());
  clouds[0].transform(transform, 0.0);
}
