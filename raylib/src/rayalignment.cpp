// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayalignment.h"
#include "simple_fft/fft.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "imagewrite.h"
#include <iostream>
#include <complex>
using namespace std;
using namespace Eigen;
using namespace RAY;

typedef complex<double> Complex;

struct Col
{
  unsigned char r,g,b,a;
};
static bool debugImageOutput = true;

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

void Array3D::operator *=(const Array3D &other)
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] *= other.cells[i];
}

void Array3D::conjugate()
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] = conj(cells[i]);
}

void Array3D::FFT()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate FFT: " << error << endl;
}

void Array3D::inverseFFT()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, dims[0], dims[1], dims[2], error))
    cout << "failed to calculate inverse FFT: " << error << endl;
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

void Array1D::operator*=(const Array1D &other)
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] *= other.cells[i];
}

void Array1D::conjugate()
{
  for (int i = 0; i<(int)cells.size(); i++)
    cells[i] = conj(cells[i]);
}

void Array1D::FFT()
{
  const char *error = nullptr;
  // TODO: change to in-place (cells not repeated)
  if (!simple_fft::FFT(*this, *this, cells.size(), error))
    cout << "failed to calculate FFT: " << error << endl;
}

void Array1D::inverseFFT()
{
  const char *error = nullptr;
  if (!simple_fft::IFFT(*this, *this, cells.size(), error))
    cout << "failed to calculate inverse FFT: " << error << endl;
}

int Array1D::maxRealIndex() const
{
  int index = 0;
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

/************************************************************************************/
void drawArray(const Array3D &array, const Vector3i &dims, const string &fileName, int index)
{
  if (!debugImageOutput)
    return;

  int width = dims[0];
  int height = dims[1];
  double maxVal = 0.0;
  for (int x = 0; x<width; x++)
  {
    for (int y = 0; y<height; y++)
    {
      double val = 0.0;
      for (int z = 0; z<dims[2]; z++)
        val += abs(array(x,y,z));
      maxVal = max(maxVal, val);
    }
  }
  
  vector<Col> pixels(width*height);
  for (int x = 0; x<width; x++)
  {
    for (int y = 0; y<height; y++)
    {
      Vector3d colour(0,0,0);
      for (int z = 0; z<dims[2]; z++)
      {
        double h = (double)z / (double)dims[2];
        Vector3d col;
        col[0] = 1.0 - h;
        col[2] = h;
        col[1] = 3.0*col[0]*col[2];
        colour += abs(array(x,y,z)) * col;
      }
      colour *= 5.0*255.0 / maxVal;
      Col col;
      col.r = clamped((int)colour[0], 0, 255);
      col.g = clamped((int)colour[1], 0, 255);
      col.b = clamped((int)colour[2], 0, 255);
      col.a = 255;
      pixels[(x+width/2)%width + width*((y+height/2)%height)] = col;
    }
  }
  stringstream str;
  str << fileName << index << ".png";
  stbi_write_png(str.str().c_str(), width, height, 4, (void *)&pixels[0], 4*width);
}

void drawArray(const vector<Array1D> &arrays, const Vector3i &dims, const string &fileName, int index)
{
  if (!debugImageOutput)
    return;

  int width = dims[0];
  int height = dims[1];
  double maxVal = 0.0;
  for (int x = 0; x<width; x++)
  {
    for (int y = 0; y<height; y++)
    {
      double val = 0.0;
      for (int z = 0; z<dims[2]; z++)
        val += abs(arrays[y + dims[1]*z](x));
      maxVal = max(maxVal, val);
    }
  }

  vector<Col> pixels(width*height);
  for (int x = 0; x<width; x++)
  {
    for (int y = 0; y<height; y++)
    {
      Vector3d colour(0,0,0);
      for (int z = 0; z<dims[2]; z++)
      {
        double h = (double)z / (double)dims[2];
        Vector3d col;
        col[0] = 1.0 - h;
        col[2] = h;
        col[1] = 3.0*col[0]*col[2];
        colour += abs(arrays[y + dims[1]*z](x)) * col;
      }
      colour *= 5.0*255.0 / maxVal;
      Col col;
      col.r = clamped((int)colour[0], 0, 255);
      col.g = clamped((int)colour[1], 0, 255);
      col.b = clamped((int)colour[2], 0, 255);
      col.a = 255;
      pixels[(x+width/2)%width + width*y] = col;
    }
  }
  stringstream str;
  str << fileName << index << ".png";
  stbi_write_png(str.str().c_str(), width, height, 4, (void *)&pixels[0], 4*width);
}

void AlignTranslationYaw::alignCloud0ToCloud1(double voxelWidth)
{
  Array3D arrays[2]; 
  // first we need to decimate the clouds into intensity grids..
  // I need to get a maximum box width, and individual boxMin, boxMaxs
  Vector3d boxMins[2], boxWidth(0,0,0);
  for (int c = 0; c<2; c++)
  {
    Vector3d boxMin(-1,-1,-1)/*1e10,1e10,1e10)*/, boxMax(-1e10,-1e10,-1e10);
    for (int i = 0; i<(int)clouds[c].ends.size(); i++)
    {
      if (clouds[c].rayBounded(i))
      {
        boxMin = minVector(boxMin, clouds[c].ends[i]);
        boxMax = maxVector(boxMax, clouds[c].ends[i]);
      }
    }
    boxMins[c] = boxMin;
    Vector3d width = boxMax - boxMin;
    boxWidth = maxVector(boxWidth, width);
  }
    
    // temporarily making the box square, to make it simple with the polar part...
 //   boxWidth = Vector3d(max(boxWidth[0], boxWidth[1]), max(boxWidth[0], boxWidth[1]), boxWidth[2]);
  
  // Now fill in the arrays with point density
  for (int c = 0; c<2; c++)
  {
    arrays[c].init(boxMins[c], boxMins[c] + boxWidth, voxelWidth);
    for (int i = 0; i<(int)clouds[c].ends.size(); i++)
      if (clouds[c].rayBounded(i))
        arrays[c](clouds[c].ends[i]) += Complex(1,0);  
    arrays[c].FFT();
    drawArray(arrays[c], arrays[c].dims, "translationInvariant", c);
  }

  Array3D &array = arrays[0];
  bool rotationToEstimate = true; // If we know there is no rotation between the clouds then we can save some cost
  if (rotationToEstimate)
  {
    // TODO: the arrays are skew symmetric so... we only need to look at half of the data....
    // however, the data is not perfectly skew symmetric for some reason...

    // OK cool, so next I need to re-map the two arrays into 4x1 grids...
    int maxRad = max(array.dims[0], array.dims[1])/2;
    Vector3i polarDims = Vector3i(4*maxRad, maxRad, array.dims[2]);
    vector<Array1D> polars[2]; // 2D grid of Array1Ds
    for (int c = 0; c<2; c++)
    {
      Array3D &a = arrays[c];
      vector<Array1D> &polar = polars[c];
      polar.resize(polarDims[1] * polarDims[2]);
      for (int j = 0; j<polarDims[1]; j++)
        for (int k = 0; k<polarDims[2]; k++)
          polar[j + polarDims[1]*k].init(polarDims[0]);

      // now map...
      for (int i = 0; i<polarDims[0]; i++)
      {
        double angle = 2.0*pi*(double)i/(double)polarDims[0];
        for (int j = 0; j<polarDims[1]; j++)
        {
          double radius = (0.5+(double)j)/(double)polarDims[1];
          Vector2d pos = radius*0.5*Vector2d((double)a.dims[0]*sin(angle), (double)a.dims[1]*cos(angle));
          if (pos[0] < 0.0)
            pos[0] += a.dims[0];
          if (pos[1] < 0.0)
            pos[1] += a.dims[1];
          int x = pos[0]; 
          int y = pos[1];
          double blendX = pos[0] - (double)x;
          double blendY = pos[1] - (double)y;
          for (int z = 0; z<polarDims[2]; z++)
          {
            // bilinear interpolation
            double val = abs(a(x,y,z)) * (1.0-blendX)*(1.0-blendY) + abs(a(x+1,y,z)) * blendX*(1.0-blendY)
                       + abs(a(x,y+1,z))*(1.0-blendX)*blendY       + abs(a(x+1,y+1,z))*blendX*blendY;
            polar[j + polarDims[1]*z](i) = Complex(radius*val, 0);
          }
        }
      }
      drawArray(polar, polarDims, "translationInvPolar", c);

      for (int j = 0; j<polarDims[1]; j++)
        for (int k = 0; k<polarDims[2]; k++)
          polar[j + polarDims[1]*k].FFT();

      drawArray(polar, polarDims, "euclideanInvariant", c);
    }

    vector<Array1D> &polar = polars[0];
    // now get the inverse fft in place:
    for (int i = 0; i<(int)polar.size(); i++)
    {
      polars[1][i].conjugate();
      polar[i] *= polars[1][i];
      polar[i].inverseFFT();
      if (i>0)
        polar[0] += polar[i]; // add all the results together into the first array
    }

    // get the angle of rotation
    int index = polar[0].maxRealIndex();
    // add a little bit of sub-pixel accuracy:
    double angle;
    int dim = polarDims[0];
    int back = (index+dim-1)%dim;
    int fwd = (index+1)%dim;
    double y0 = polar[0](back).real();
    double y1 = polar[0](index).real();
    double y2 = polar[0](fwd).real();
    angle = index + 0.5*(y0 - y2)/(y0 + y2 - 2.0*y1); // just a quadratic maximum -b/2a for heights y0,y1,y2
    // but the FFT wraps around, so:
    if (angle > dim / 2)
      angle -= dim;
    angle *= 2.0*pi/(double)polarDims[0];
    cout << "estimated yaw rotation: " << angle << endl;

    // ok, so let's rotate A towards B, and re-run the translation FFT
    clouds[0].transform(Pose(Vector3d(0,0,0), Quaterniond(AngleAxisd(angle, Vector3d(0,0,1)))), 0.0);

    boxMins[0] = Vector3d(1e10,1e10,1e10);
    for (int i = 0; i<(int)clouds[0].ends.size(); i++)
      if (clouds[0].rayBounded(i))
        boxMins[0] = minVector(boxMins[0], clouds[0].ends[i]);
    arrays[0].cells.clear();
    arrays[0].init(boxMins[0], boxMins[0]+boxWidth, voxelWidth);

    for (int i = 0; i<(int)clouds[0].ends.size(); i++)
      if (clouds[0].rayBounded(i))
        arrays[0](clouds[0].ends[i]) += Complex(1,0);  // TODO: this could go out of bounds!
    arrays[0].FFT();
  }

  /****************************************************************************************************/
  // now get the the translation part
  arrays[1].conjugate();
  arrays[0] *= arrays[1];
  arrays[0].inverseFFT();

  // find the peak
  Vector3i ind = array.maxRealIndex();
  // add a little bit of sub-pixel accuracy:
  Vector3d pos;
  for (int axis = 0; axis<3; axis++)
  {
    Vector3i back = ind, fwd = ind;
    int &dim = array.dims[axis];
    back[axis] = (ind[axis]+dim-1)%dim;
    fwd[axis] = (ind[axis]+1)%dim;
    double y0 = array(back).real();
    double y1 = array(ind).real();
    double y2 = array(fwd).real();
    pos[axis] = ind[axis] + 0.5*(y0 - y2)/(y0 + y2 - 2.0*y1); // just a quadratic maximum -b/2a for heights y0,y1,y2
    // but the FFT wraps around, so:
    if (pos[axis] > dim / 2)
      pos[axis] -= dim;
  }
  pos *= -array.voxelWidth;
  pos += boxMins[1]-boxMins[0];
  cout << "estimated translation: " << pos.transpose() << endl;

  Pose transform(pos, Quaterniond::Identity());
  clouds[0].transform(transform, 0.0);
}