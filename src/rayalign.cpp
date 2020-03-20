// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raydraw.h"
#include "rayalignment.h"
#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <complex>
using namespace std;
using namespace Eigen;
using namespace RAY;
typedef complex<double> Complex;


void usage(bool error=false)
{
  cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << endl;
  cout << "usage:" << endl;
  cout << "rayalign raycloudA raycloudB." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "rayalign");
  DebugDraw draw;
  if (argc != 3)
    usage();

  const double voxelWidth = 0.5;
  string fileA = argv[1];
  string fileB = argv[2];
  Cloud clouds[2];
  clouds[0].load(fileA);
  clouds[1].load(fileB);

  Array3D arrays[2]; 
  // first we need to decimate the clouds into intensity grids..
  // I need to get a maximum box width, and individual boxMin, box Maxs
  Vector3d boxMins[2], boxWidth(0,0,0);
  for (int c = 0; c<2; c++)
  {
    Vector3d boxMin(-1,-1,-1)/*1e10,1e10,1e10)*/, boxMax(-1e10,-1e10,-1e10);
    for (auto &point: clouds[c].ends)
    {
      boxMin = minVector(boxMin, point);
      boxMax = maxVector(boxMax, point);
    }
    boxMins[c] = boxMin;
    Vector3d width = boxMax - boxMin;
    boxWidth = maxVector(boxWidth, width);
  }
    
    // temporarily making the box square, to make it simple with the polar part...
    boxWidth = Vector3d(max(boxWidth[0], boxWidth[1]), max(boxWidth[0], boxWidth[1]), boxWidth[2]);
  
  // Now fill in the arrays with point density
  for (int c = 0; c<2; c++)
  {
    arrays[c].init(boxMins[c], boxMins[c] + boxWidth, voxelWidth);
    for (auto &point: clouds[c].ends)
      arrays[c](point) += Complex(1,0);  
    arrays[c].fft();
  }
  Array3D &array = arrays[0];

  bool rotationToEstimate = true; // If we know there is no rotation between the clouds then we can save some cost
  if (rotationToEstimate)
  {
    // OK cool, so next I need to re-map the two arrays into 2x1 grids...
    int maxRad = max(array.dims[0], array.dims[1])/2;
    Vector3i polarDims = Vector3i(4*maxRad, maxRad, array.dims[2]);
    vector<Array1D> polars[2]; // 2D grid of Array1Ds
    for (int c = 0; c<2; c++)
    {
      Array3D &a = arrays[c];
      vector<Array1D> &polar = polars[c];
      polar.resize(maxRad * array.dims[2]);
      for (int j = 0; j<polarDims[1]; j++)
        for (int k = 0; k<polarDims[2]; k++)
          polar[j + polarDims[1]*k].init(polarDims[0]);

      // now map...
      for (int i = 0; i<polarDims[0]; i++)
      {
        double angle = 2.0*pi*(double)i/(double)polarDims[0];
        for (int j = 0; j<polarDims[1]; j++)
        {
          double radius = (double)maxRad * (double)j/(double)polarDims[1];
          Vector2d pos = radius * Vector2d(sin(angle), cos(angle));
          if (pos[0] < 0.0)
            pos[0] += 2.0*maxRad;
          if (pos[1] < 0.0)
            pos[1] += 2.0*maxRad;
          int x = pos[0]; 
          int y = pos[1];
          double blendX = pos[0] - (double)x;
          double blendY = pos[1] - (double)y;
          for (int z = 0; z<polarDims[2]; z++)
          {
            // bilinear interpolation
            double val = a(x,y,z).real() * (1.0-blendX)*(1.0-blendY) + a(x+1,y,z).real() * blendX*(1.0-blendY)
                      + a(x,y+1,z).real()*(1.0-blendX)*blendY       + a(x+1,y+1,z).real()*blendX*blendY;
            polar[j + polarDims[1]*z](i) = Complex(radius*val, 0);
          }
        }
      }
      for (int j = 0; j<polarDims[1]; j++)
        for (int k = 0; k<polarDims[2]; k++)
          polar[j + polarDims[1]*k].fft();
    }

    vector<Array1D> &polar = polars[0];
    // now get the inverse fft in place:
    for (int i = 0; i<(int)polar.size(); i++)
    {
      polar[i].cwiseProductInplace(polars[1][i].conjugateInplace()).ifft(); 
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
    angle *= -2.0*pi/(double)polarDims[0];
    cout << "found angle: " << angle << endl;

    // ok, so let's rotate A towards B, and re-run the translation FFT
    clouds[0].transform(Pose(Vector3d(0,0,0), Quaterniond(AngleAxisd(angle, Vector3d(0,0,1)))), 0.0);

    boxMins[0] = Vector3d(1e10,1e10,1e10);
    for (auto &point: clouds[0].ends)
      boxMins[0] = minVector(boxMins[0], point);
    arrays[0].cells.clear();
    arrays[0].init(boxMins[0], boxMins[0]+boxWidth, voxelWidth);

    for (auto &point: clouds[0].ends)
      arrays[0](point) += Complex(1,0);  // TODO: this could go out of bounds!
    arrays[0].fft();
  }

  /****************************************************************************************************/
  // now get the the translation part
  arrays[0].cwiseProductInplace(arrays[1].conjugateInplace()).ifft();  

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
  cout << "translation: " << pos.transpose() << " plus boxMin difference: " << (boxMins[1]-boxMins[0]).transpose() << " gives: " << (pos + boxMins[1]-boxMins[0]).transpose() << endl;
  pos += boxMins[1]-boxMins[0];

  Pose transform(pos, Quaterniond::Identity());
  clouds[0].transform(transform, 0.0);

  string fileStub = fileA;
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);
  clouds[0].save(fileStub + "_aligned.ply");  

  return true;
}
