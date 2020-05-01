// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace ray;

void usage(int exit_code = 0)
{
  cout << "Colour the ray cloud, and/or shade it" << endl;
  cout << "usage:" << endl;
  cout << "raycolour raycloud time shaded  - colour by time, and shaded (optional on all types)" << endl;
  cout << "                   height       - colour by height" << endl;
  cout << "                   shape        - colour by geometry shape (r,g,b: spherical, cylinderical, planar)" << endl;
  cout << "                   normal       - colour by normal" << endl;
  cout << "                   1,1,1 shaded - just (r,g,b) shaded" << endl;
  exit(exit_code);
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
  const int search_size = 20;
  vector<Vector3d> centroids;
  vector<Vector3d> dimensions;
  vector<Vector3d> normals; 
  vector<Matrix3d> matrices;
  MatrixXi indices;
  vector<Vector3d> *cents = NULL, *dims = NULL, *norms = NULL;
  vector<Matrix3d> *mats = NULL;
  MatrixXi *inds = NULL;

  // what do we want to calculate...
  bool calc_surfels = true;
  if (type == "normal")
    norms = &normals;
  else if (type == "shape")
    dims = &dimensions;
  else
    calc_surfels = shaded;
  if (shaded)
  {
    norms = &normals;
    inds = &indices;
    cents = &centroids;
  }

  if (calc_surfels)
    cloud.getSurfels(search_size, cents, norms, dims, mats, inds);
  
  // Q: can I do better? in particular, can I colour as doubles and only quantise at the end?
  if (type.find(",") != string::npos)
  {
    stringstream ss(type);
    Vector3d vec;
    ss >> vec[0]; ss.ignore(1); ss>>vec[1]; ss.ignore(1); ss>>vec[2];
    for (auto &col: cloud.colours)
    {
      col.red = (uint8_t)(255.0*vec[0]);
      col.green = (uint8_t)(255.0*vec[1]);
      col.blue = (uint8_t)(255.0*vec[2]);
    }
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
      cloud.colours[i].red   = (uint8_t)(255.0*(0.3 + 0.7*sphericity));
      cloud.colours[i].green = (uint8_t)(255.0*(0.3 + 0.7*cylindricality)); 
      cloud.colours[i].blue  = (uint8_t)(255.0*(0.3 + 0.7*planarity));
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
      double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0, sum_yy = 0, n = 0;
      for (int j = 0; j<search_size && indices(j,i)>-1; j++)
      {
        int id = indices(j,i);
        Vector3d flat = cloud.ends[id] - centroids[i];
        double y = flat.dot(normals[i]);
        flat -= y*normals[i];
        double x = flat.squaredNorm();
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_xx += x*x;
        sum_yy += y*y;
        n++;
      }
      double den = n*sum_xx - sum_x*sum_x;
      if (abs(den) < 1e-8)
        curvatures[i] = 0.0;
      else 
        curvatures[i] = (n*sum_xy - sum_x*sum_y) / den;
    }
    Vector3d light_dir = Vector3d(0.2, 0.4, 1.0).normalized();
    double curve_scale = 4.0;
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double scale1 = 0.5 + 0.5*normals[i].dot(light_dir);
      double scale2 = 0.5 - 0.5*curvatures[i]/curve_scale; 
      double s = 0.25 + 0.75*clamped((scale1 + scale2)/2.0, 0.0, 1.0);
      cloud.colours[i].red = (uint8_t)((double)cloud.colours[i].red * s);
      cloud.colours[i].green = (uint8_t)((double)cloud.colours[i].green * s);
      cloud.colours[i].blue = (uint8_t)((double)cloud.colours[i].blue * s);
    }
  }

  string file_stub = file;
  if (file.substr(file.length()-4)==".ply")
    file_stub = file.substr(0,file.length()-4);
  cloud.save(file_stub + "_coloured.ply");
  return true;
}
