// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Colour the ray cloud, and/or shade it" << endl;
  cout << "usage:" << endl;
  cout << "raycolour raycloud time shaded  - colour by time, and shaded (optional on all types)" << endl;
  cout << "                   height       - colour by height" << endl;
  cout << "                   shape        - colour by geometry shape (r,g,b: spherical, cylinderical, planar)" << endl;
  cout << "                   normal       - colour by normal" << endl;
  cout << "                   white shaded - just shaded" << endl;
  exit(error);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  if (argc != 3 && argc != 4)
    usage();
  
  bool shaded = false;
  if (argc == 4)
  {
    shaded = true;
    if (string(argv[3]) != "shaded")
      usage();
  }
  string file = argv[1];
  Cloud cloud;
  cloud.load(file);
  string type = string(argv[2]);

  // what I need is the normal, curvature, eigenvalues, per point. 
  struct Data
  {
    Vector3d normal;
    Vector3d dimensions;
    double curvature;
  };
  vector<Data> data(cloud.ends.size());
  const int searchSize = 20;
  vector<Vector3d> centroids;
  vector<Vector3d> dimensions;
  vector<Vector3d> normals; 
  vector<Matrix3d> matrices;
  MatrixXi indices;
  vector<Vector3d> *cents = NULL, *dims = NULL, *norms = NULL;
  vector<Matrix3d> *mats = NULL;
  MatrixXi *inds = NULL;

  // what do we want to calculate...
  bool calcSurfels = true;
  if (type == "normal")
    norms = &normals;
  else if (type == "shape")
    dims = &dimensions;
  else
    calcSurfels = shaded;
  if (shaded)
  {
    norms = &normals;
    inds = &indices;
    cents = &centroids;
  }

  if (calcSurfels)
    cloud.getSurfels(searchSize, cents, norms, dims, mats, inds);
  
  // Q: can I do better? in particular, can I colour as doubles and only quantise at the end?
  if (type == "white")
  {
    for (auto &col: cloud.colours)
      col.red = col.green = col.blue = 255;
  }
  else if (type == "time")
  {
    colourByTime(cloud.times, cloud.colours, false);
  }
  else if (type == "height")
  {
    vector<double> heights(cloud.ends.size());
    for (int i = 0; i<(int)cloud.ends.size(); i++)
      heights[i] = cloud.ends[i][2];
    redGreenBlueSpectrum(heights, cloud.colours, 10.0, false);
  }
  else if (type == "shape")
  {
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double sphericity = dimensions[i][0] / dimensions[i][2];
      double cylindricality = 1.0 - dimensions[i][1] / dimensions[i][2];
      double planarity = 1.0 - dimensions[i][0] / dimensions[i][1];
      cloud.colours[i].red = (uint8_t)(255.0*sphericity);
      cloud.colours[i].green = (uint8_t)(255.0*cylindricality);
      cloud.colours[i].blue = (uint8_t)(255.0*planarity);
    }
  }
  else if (type == "normal")
  {
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      cloud.colours[i].red = (uint8_t)(255.0*(0.5 + 0.5*normals[i][0]));
      cloud.colours[i].green = (uint8_t)(255.0*(0.5 + 0.5*normals[i][1]));
      cloud.colours[i].blue = (uint8_t)(255.0*(0.5 + 0.5*normals[i][2]));
    }
  }
  else
    usage();

  if (shaded)
  {
    vector<double> curvatures(cloud.ends.size());
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0, n = 0;
      for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
      {
        int id = indices(j,i);
        Vector3d flat = cloud.ends[id] - centroids[i];
        double y = flat.dot(normals[i]);
        flat -= y*normals[i];
        double x = flat.squaredNorm();
        sumX += x;
        sumY += y;
        sumXY += x*y;
        sumXX += x*x;
        sumYY += y*y;
        n++;
      }
      double den = n*sumXX - sumX*sumX;
      if (abs(den) < 1e-8)
        curvatures[i] = 0.0;
      else 
        curvatures[i] = (n*sumXY - sumX*sumY) / den;
    }
    Vector3d lightDir = Vector3d(0.2, 0.4, 1.0).normalized();
    double curveScale = 4.0;
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double scale1 = 0.5 + 0.5*normals[i].dot(lightDir);
      double scale2 = 0.5 + 0.5*clamped(-curvatures[i]/curveScale, -1.0, 1.0); 
      double s = 0.25 + 0.75*(scale1 + scale2)/2.0;
      cloud.colours[i].red = (uint8_t)((double)cloud.colours[i].red * s);
      cloud.colours[i].green = (uint8_t)((double)cloud.colours[i].green * s);
      cloud.colours[i].blue = (uint8_t)((double)cloud.colours[i].blue * s);
    }
  }

  string fileStub = file;
  if (file.substr(file.length()-4)==".ply")
    fileStub = file.substr(0,file.length()-4);
  cloud.save(fileStub + "_coloured.ply");
  return true;
}
