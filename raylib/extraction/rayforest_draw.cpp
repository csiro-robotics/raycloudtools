// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayforest.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imagewrite.h"

namespace ray
{
struct Col
{
  Col(){}
  Col(uint8_t shade) : r(shade), g(shade), b(shade), a(255) {}
  void operator +=(const Col &col)
  {
    r = (uint8_t)std::min((int)r + (int)col.r, 255);
    g = (uint8_t)std::min((int)g + (int)col.g, 255);
    b = (uint8_t)std::min((int)b + (int)col.b, 255);
    a = (uint8_t)std::min((int)a + (int)col.a, 255);
  }
  uint8_t r, g, b, a;
};

void Forest::drawSegmentation(const std::string &filename, const std::vector<TreeNode> &trees)
{
  const double max_tree_height = 30.0; // TODO: fix!
  if (!verbose)
    return;
  Field2D<Col> pixels((int)indexfield_.rows(), (int)indexfield_.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
  {
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      int ind = indexfield_(x, y);
      if (ind == -1)
        pixels(x, y) = Col(0);
      else
      {
        while (trees[ind].attaches_to != -1)
          ind = trees[ind].attaches_to;
        srand(1 + ind);
        Col col;
        col.a = 255;
        col.r = (uint8_t)(rand()%256);
        col.g = (uint8_t)(rand()%256);
        col.b = (uint8_t)(rand()%256);
        pixels(x, y) = col;
      }
    }
  }
 /* for (int i = 0; i<(int)trees.size(); i++)
  {
    int ind = i;
    while (trees[ind].attaches_to != -1)
      ind = trees[ind].attaches_to;
    const TreeNode &tree = trees[ind];
    Eigen::Vector2i max_bound = tree.max_bound;
    Eigen::Vector2i d = tree.max_bound - tree.min_bound;
    if (d[0] < 0)
      max_bound[0] += res;
    if (d[1] < 0)
      max_bound[1] += res;
    for (int x = tree.min_bound[0]; x <= max_bound[0]; x++)
    {
      pixels[(x%res) + res * tree.min_bound[1]] = Col(255);
      pixels[(x%res) + res * (max_bound[1]%res)] = Col(255);
    }
    for (int y = tree.min_bound[1]; y <= max_bound[1]; y++)
    {
      pixels[tree.min_bound[0] + res * (y%res)] = Col(255);
      pixels[(max_bound[0]%res) + res * (y%res)] = Col(255);
    }
  }*/
  for (int i = 0; i<(int)trees.size(); i++)
  {
    int ind = i;
    while (trees[ind].attaches_to != -1)
      ind = trees[ind].attaches_to;
    const TreeNode &tree = trees[ind];

    Vector4d abcd = tree.curv_mat.ldlt().solve(tree.curv_vec);
    double a = abcd[0], b = abcd[1], c = abcd[2], d = abcd[3];
    double x = -b/(2*a);
    double y = -c/(2*a);
    double z = d - (b*b + c*c)/(4*a);
    int X = (int)x; 
    int Y = (int)y; 
    if (X>=0 && X<pixels.dims[0] && Y >= 0.0 && Y<pixels.dims[1])
      pixels(X, Y) = Col((uint8_t)(255.0 * ray::clamped(z/max_tree_height, 0.0, 1.0)));
  }

  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield)
{
  if (!verbose)
    return;
  double max_height = -1e10; // heightfield.maxCoeff();
  for (int i = 0; i<heightfield.rows(); i++)
    for (int j = 0; j<heightfield.cols(); j++)
      if (heightfield(i,j) < 1000.0)
        max_height = std::max(max_height, heightfield(i,j));

  Field2D<Col> pixels((int)heightfield.rows(), (int)heightfield.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
    for (int y = 0; y < pixels.dims[1]; y++)
      pixels(x, y) = Col((uint8_t)(255.0 * heightfield(x, y)/max_height));
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawGraph(const std::string &filename, const std::vector<Eigen::Vector3d> &data, double x_max, double y_max, double strength_max)
{
  if (!verbose)
    return;
  int res = 256;
  Field2D<Col> pixels(res, res);
  for (auto &c: pixels.data)
    c = Col(20);
  for (auto &item: data)
  {
    double x = (double)(res - 1) * item[0] / x_max;
    double y = (double)(res - 1) * item[1] / y_max;
    if (x >= 0.0 && x<(double)res-1.0 && y >= 0.0 && y<(double)res-1.0)
      pixels((int)x, (int)y) = Col((uint8_t)(255.0 * item[2]/strength_max));
  }
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawTrees(const std::string &filename, const Forest::Result &result)
{
  const double max_height = 30.0; // TODO: fix!
  // I should probably draw the result
  if (!verbose)
    return;
  const int res = 256; // need to get actual bounds
  Field2D<Col> pixels(res, res);
  for (auto &c: pixels.data)
    c = Col(0); 
  for (auto &tip: result.tree_tips)
  {
    double length = tip[2] - result.ground_height;
    double crown_radius = length/result.treelength_per_crownradius;
    double curvature = 1.0 / crown_radius;
    double draw_radius = std::min(1.2 * crown_radius, 50.0); 
    for (int x = (int)(tip[0] - draw_radius); x<= (int)(tip[0]+draw_radius); x++)
    {
      for (int y = (int)(tip[1] - draw_radius); y<= (int)(tip[1]+draw_radius); y++)
      {
        if (x < 0 || x >= res || y<0 || y>=res)
          continue;
        double X = (double)x - tip[0];
        double Y = (double)y - tip[1];
        double mag2 = (double)(X*X + Y*Y);
        if (mag2 <= draw_radius*draw_radius)
        {
          double height = tip[2] - mag2 * curvature;
          double shade = (height - result.ground_height)/(max_height - result.ground_height);
          Col col(uint8_t(255.0*shade));
          if (pixels(x, y).r < col.r)
            pixels(x, y) = col;
        }
      }
    }
  }    
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}
} // namespace ray