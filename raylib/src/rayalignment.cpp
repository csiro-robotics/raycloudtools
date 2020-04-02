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
// #define WEIGHTED_ALIGN

struct Col
{
  unsigned char r,g,b,a;
};

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
  nullCell = 0;
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

void Array3D::fillWithRays(const Cloud &cloud)
{
  // unlike the end point densities, the weight is just 0 or 1, but requires walking through the grid for every ray
  // maybe a better choice would be a reuseable 'volume' function (occupancy grid).
  for (int i = 0; i<(int)cloud.ends.size(); i++)
  {
    // TODO: Should be a lmbda function!
    Vector3d dir = cloud.ends[i] - cloud.starts[i];
    Vector3d dirSign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (cloud.starts[i] - boxMin)/voxelWidth;
    Vector3d end = (cloud.ends[i] - boxMin)/voxelWidth;
    Vector3i startIndex(start[0], start[1], start[2]);
    Vector3i endIndex(end[0], end[1], end[2]);
    double lengthSqr = (endIndex - startIndex).squaredNorm();
    Vector3i index = startIndex;
    while ((index - startIndex).squaredNorm()<=lengthSqr+1e-10)
    {
      if (index[0] >= 0 && index[0]<dims[0] && index[1] >= 0 && index[1]<dims[1] && index[2] >= 0 && index[2]<dims[2])
        (*this)(index[0], index[1], index[2]) += Complex(1,0); // add weight to these areas...

      Vector3d mid = boxMin + voxelWidth*Vector3d(index[0]+0.5, index[1]+0.5, index[2]+0.5);
      Vector3d nextBoundary = mid + 0.5*voxelWidth*dirSign;
      Vector3d delta = nextBoundary - cloud.starts[i];
      Vector3d d(delta[0]/dir[0], delta[1]/dir[1], delta[2]/dir[2]);
      if (d[0] < d[1] && d[0]<d[2])
        index[0] += dirSign[0];
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += dirSign[1];
      else
        index[2] += dirSign[2];
    }
  }
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

void drawArray(const Array3D &array, const Vector3i &dims, const string &fileName, int index)
{
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

void Array1D::polarCrossCorrelation(const Array3D *arrays, bool verbose, double *maxWeights)
{
  // OK cool, so next I need to re-map the two arrays into 4x1 grids...
  int maxRad = max(arrays[0].dims[0], arrays[0].dims[1])/2;
  Vector3i polarDims = Vector3i(4*maxRad, maxRad, arrays[0].dims[2]);
  vector<Array1D> polars[2];
  for (int c = 0; c<2; c++)
  {
    const Array3D &a = arrays[c];
    vector<Array1D> &polar = polars[c];
    polar.resize(polarDims[1] * polarDims[2]);
    for (int j = 0; j<polarDims[1]; j++)
      for (int k = 0; k<polarDims[2]; k++)
        polar[j + polarDims[1]*k].init(polarDims[0]);

    // now map...
    for (int i = 0; i<polarDims[0]; i++)
    {
      double angle = 2.0*pi*(double)(i+0.5)/(double)polarDims[0];
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
        int x2 = (x+1)%a.dims[0];
        int y2 = (y+1)%a.dims[1];
        double blendX = pos[0] - (double)x;
        double blendY = pos[1] - (double)y;
        for (int z = 0; z<polarDims[2]; z++)
        {
          // bilinear interpolation
          double val = abs(a(x,y,z)) * (1.0-blendX)*(1.0-blendY) + abs(a(x2,y,z)) * blendX*(1.0-blendY)
                      + abs(a(x,y2,z))*(1.0-blendX)*blendY       + abs(a(x2,y2,z))*blendX*blendY;
          polar[j + polarDims[1]*z](i) = Complex(radius*val, 0);
        }
      }
    }

    if (maxWeights)
    {
      maxWeights[c] = 0;
      for (int j = 0; j<polarDims[1]; j++)
        for (int k = 0; k<polarDims[2]; k++)
          for (int l = 0; l<(int)polar[j + polarDims[1]*k].cells.size(); l++)
            maxWeights[c] += norm(polar[j + polarDims[1]*k].cells[l]);
    }
    if (verbose)
      drawArray(polar, polarDims, "translationInvPolar", c);
    for (int j = 0; j<polarDims[1]; j++)
      for (int k = 0; k<polarDims[2]; k++)
        polar[j + polarDims[1]*k].FFT();
    if (verbose)
      drawArray(polar, polarDims, "euclideanInvariant", c);
  }

  // now get the inverse fft in place:
  init(polarDims[0]);
  for (int i = 0; i<(int)polars[0].size(); i++)
  {
    polars[1][i].conjugate();
    polars[0][i] *= polars[1][i];
    polars[0][i].inverseFFT();
    (*this) += polars[0][i]; // add all the results together into the first array
  }
}

/************************************************************************************/

void AlignTranslationYaw::alignCloud0ToCloud1(double voxelWidth, bool verbose)
{
  // first we need to decimate the clouds into intensity grids..
  // I need to get a maximum box width, and individual boxMin, boxMaxs
  Vector3d boxMins[2], boxWidth(0,0,0);
  for (int c = 0; c<2; c++)
  {
    Vector3d boxMin(1e10,1e10,1e10), boxMax(-1e10,-1e10,-1e10);
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
  
  Array3D arrays[2]; 
#if defined(WEIGHTED_ALIGN)
  // very small below 1 and this will find small areas of overlap, very large above 1 and it is
  // like the unweighted Fourier-Mellin, which maximises the overlap for all points
  const double stability = 10.0; 
  double maxWeights[2] = {0,0};
  Array3D weights[2];
#endif
  // Now fill in the arrays with point density
  for (int c = 0; c<2; c++)
  {
    arrays[c].init(boxMins[c], boxMins[c] + boxWidth, voxelWidth);
    for (int i = 0; i<(int)clouds[c].ends.size(); i++)
      if (clouds[c].rayBounded(i))
        arrays[c](clouds[c].ends[i]) += Complex(1,0);  
#if defined(WEIGHTED_ALIGN)
    weights[c].init(boxMins[c], boxMins[c] + boxWidth, voxelWidth);
    weights[c].fillWithRays(clouds[c]);
    for (auto &comp: weights[c].cells)
      maxWeights[c] += norm(comp);
    arrays[c] *= weights[c];
    weights[c].FFT();
    if (verbose) 
      drawArray(weights[c], weights[c].dims, "translationInvariant_weight", c);
#endif
    arrays[c].FFT();
    if (verbose) 
      drawArray(arrays[c], arrays[c].dims, "translationInvariant", c);
  }

  bool rotationToEstimate = true; // If we know there is no rotation between the clouds then we can save some cost
  if (rotationToEstimate)
  {
    Array1D polar; 
    polar.polarCrossCorrelation(arrays, verbose);
    #if defined(WEIGHTED_ALIGN)
    Array1D polarWeight; 
    double maxPolarWeights[2] = {0,0};
    polarWeight.polarCrossCorrelation(weights, verbose, maxPolarWeights);
    polar.stableDivideBy(polarWeight, stability*min(maxPolarWeights[0], maxPolarWeights[1]));
    #endif

    // get the angle of rotation
    int index = polar.maxRealIndex();
    // add a little bit of sub-pixel accuracy:
    double angle;
    int dim = polar.cells.size();
    int back = (index+dim-1)%dim;
    int fwd = (index+1)%dim;
    double y0 = polar(back).real();
    double y1 = polar(index).real();
    double y2 = polar(fwd).real();
    angle = index + 0.5*(y0 - y2)/(y0 + y2 - 2.0*y1); // just a quadratic maximum -b/2a for heights y0,y1,y2
    // but the FFT wraps around, so:
    if (angle > dim / 2)
      angle -= dim;
    angle *= 2.0*pi/(double)polar.cells.size();
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
        arrays[0](clouds[0].ends[i]) += Complex(1,0);  
    #if defined(WEIGHTED_ALIGN)
    weights[0].cells.clear();
    weights[0].init(boxMins[0], boxMins[0]+boxWidth, voxelWidth);    
    weights[0].fillWithRays(clouds[0]);
    maxWeights[0] = 0;
    for (auto &comp: weights[0].cells)
      maxWeights[0] += norm(comp);
    arrays[0] *= weights[0];
    weights[0].FFT();
    #endif
    arrays[0].FFT();
  }

  /****************************************************************************************************/
  // now get the the translation part
  arrays[1].conjugate();
  arrays[0] *= arrays[1];
  arrays[0].inverseFFT();
  #if defined(WEIGHTED_ALIGN)
  weights[1].conjugate();
  weights[0] *= weights[1];
  weights[0].inverseFFT();  
  arrays[0].stableDivideBy(weights[0], stability*min(maxWeights[0], maxWeights[1])); // TODO: what should this factor be??
  #endif

  // find the peak
  Array3D &array = arrays[0];
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