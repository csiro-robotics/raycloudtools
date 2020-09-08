// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayparse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage(int exit_code = 0)
{
  std::cout << "Colour the ray cloud, and/or shade it" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycolour raycloud time          - colour by time (optional on all types)" << std::endl;
  std::cout << "                   height        - colour by height" << std::endl;
  std::cout << "                   shape         - colour by geometry shape (r,g,b: spherical, cylinderical, planar)" << std::endl;
  std::cout << "                   normal        - colour by normal" << std::endl;
  std::cout << "                   alpha         - colour by alpha channel (which typically represents intensity)" << std::endl;
  std::cout << "                   1,1,1         - just (r,g,b)" << std::endl;
  std::cout << "                         --unlit - flat shaded" << std::endl;
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  ray::KeyChoice colour_type({"time", "height", "shape", "normal", "alpha"});
  ray::OptionalFlagArgument unlit("unlit", 'u');
  ray::Vector3dArgument col;
  bool flat_colour = ray::parseCommandLine(argc, argv, {&cloud_file, &col, &unlit});
  if (!flat_colour && !ray::parseCommandLine(argc, argv, {&cloud_file, &colour_type, &unlit}))
    usage();
  
  bool shading = !unlit.is_set;
  ray::Cloud cloud;
  if (!cloud.load(cloud_file.name))
    return 0;
  std::string type = colour_type.selected_key;

  // what I need is the normal, curvature, eigenvalues, per point.
  struct Data
  {
    Eigen::Vector3d normal;
    Eigen::Vector3d dimensions;
    double curvature;
  };
  std::vector<Data> data(cloud.ends.size());
  const int search_size = 20;
  std::vector<Eigen::Vector3d> centroids;
  std::vector<Eigen::Vector3d> dimensions;
  std::vector<Eigen::Vector3d> normals;
  Eigen::MatrixXi indices;
  std::vector<Eigen::Vector3d> *cents = NULL, *dims = NULL, *norms = NULL;
  std::vector<Eigen::Matrix3d> *mats = NULL;
  Eigen::MatrixXi *inds = NULL;

  // what do we want to calculate...
  bool calc_surfels = true;
  if (type == "normal")
    norms = &normals;
  else if (type == "shape")
    dims = &dimensions;
  else
    calc_surfels = shading;
  if (shading)
  {
    norms = &normals;
    inds = &indices;
    cents = &centroids;
  }

  if (calc_surfels)
    cloud.getSurfels(search_size, cents, norms, dims, mats, inds);

  if (flat_colour)
  {
    for (auto &colour : cloud.colours)
    {
      colour.red = (uint8_t)(255.0 * col.value[0]);
      colour.green = (uint8_t)(255.0 * col.value[1]);
      colour.blue = (uint8_t)(255.0 * col.value[2]);
    }
  }
  else if (type == "time")
  {
    colourByTime(cloud.times, cloud.colours, false);
  }
  else if (type == "height")
  {
    std::vector<double> heights(cloud.ends.size());
    for (int i = 0; i < (int)cloud.ends.size(); i++) heights[i] = cloud.ends[i][2];
    redGreenBlueSpectrum(heights, cloud.colours, 10.0, false);
  }
  else if (type == "shape")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double sphericity = dimensions[i][0] / dimensions[i][2];
      double cylindricality = 1.0 - dimensions[i][1] / dimensions[i][2];
      double planarity = 1.0 - dimensions[i][0] / dimensions[i][1];
      cloud.colours[i].red = (uint8_t)(255.0 * (0.3 + 0.7 * sphericity));
      cloud.colours[i].green = (uint8_t)(255.0 * (0.3 + 0.7 * cylindricality));
      cloud.colours[i].blue = (uint8_t)(255.0 * (0.3 + 0.7 * planarity));
    }
  }
  else if (type == "normal")
  {
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      cloud.colours[i].red = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][0]));
      cloud.colours[i].green = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][1]));
      cloud.colours[i].blue = (uint8_t)(255.0 * (0.5 + 0.5 * normals[i][2]));
    }
  }
  else if (type == "alpha")
  {
    std::vector<double> alphas(cloud.colours.size());
    for (int i = 0; i < (int)cloud.colours.size(); i++)
      alphas[i] = double(cloud.colours[i].alpha);
    const bool replace_alpha = false;
    redGreenBlueGradient(alphas, cloud.colours, 0.0, 255.0, replace_alpha);
  }
  else
    usage();

  if (shading)
  {
    std::vector<double> curvatures(cloud.ends.size());
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0, sum_yy = 0, n = 0;
      for (int j = 0; j < search_size && indices(j, i) > -1; j++)
      {
        int id = indices(j, i);
        Eigen::Vector3d flat = cloud.ends[id] - centroids[i];
        double y = flat.dot(normals[i]);
        flat -= y * normals[i];
        double x = flat.squaredNorm();
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        n++;
      }
      double den = n * sum_xx - sum_x * sum_x;
      if (abs(den) < 1e-8)
        curvatures[i] = 0.0;
      else
        curvatures[i] = (n * sum_xy - sum_x * sum_y) / den;
    }
    Eigen::Vector3d light_dir = Eigen::Vector3d(0.2, 0.4, 1.0).normalized();
    double curve_scale = 4.0;
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (!cloud.rayBounded(i))
        continue;
      double scale1 = 0.5 + 0.5 * normals[i].dot(light_dir);
      double scale2 = 0.5 - 0.5 * curvatures[i] / curve_scale;
      double s = 0.25 + 0.75 * ray::clamped((scale1 + scale2) / 2.0, 0.0, 1.0);
      cloud.colours[i].red = (uint8_t)((double)cloud.colours[i].red * s);
      cloud.colours[i].green = (uint8_t)((double)cloud.colours[i].green * s);
      cloud.colours[i].blue = (uint8_t)((double)cloud.colours[i].blue * s);
    }
  }

  cloud.save(cloud_file.nameStub() + "_coloured.ply");
  return true;
}
