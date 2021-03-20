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
  Col(uint8_t red, uint8_t green, uint8_t blue, uint8_t alpha) : r(red), g(green), b(blue), a(alpha) {}
  void operator +=(const Col &col)
  {
    r = (uint8_t)std::min((int)r + (int)col.r, 255);
    g = (uint8_t)std::min((int)g + (int)col.g, 255);
    b = (uint8_t)std::min((int)b + (int)col.b, 255);
    a = (uint8_t)std::min((int)a + (int)col.a, 255);
  }
  uint8_t r, g, b, a;
};

void Forest::drawFinalSegmentation(const std::string &filename, std::vector<TreeNode> &trees, std::vector<int> &indices)
{
  if (!verbose)
    return;
  Field2D<Col> pixels((int)indexfield_.rows(), (int)indexfield_.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
  {
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      int ind = indexfield_(x, y);
      Col col;
      if (ind == -1)
        pixels(x, y) = Col(0);
      else
      {
        while (trees[ind].attaches_to != -1)
        {
          if (std::find(indices.begin(), indices.end(), ind) != indices.end())
            break;
          ind = trees[ind].attaches_to;
        }
        if (std::find(indices.begin(), indices.end(), ind) == indices.end())
          continue;
        srand(1 + ind);
        col.a = 255;
        col.r = (uint8_t)(rand()%256);
        col.g = (uint8_t)(rand()%256);
        col.b = (uint8_t)(rand()%256);
        pixels(x, y) = col;
      }
    }
  }
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}


void Forest::drawSegmentation(const std::string &filename, std::vector<TreeNode> &trees)
{
  if (!verbose)
    return;
  Field2D<Col> pixels((int)indexfield_.rows(), (int)indexfield_.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
  {
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      Eigen::Vector3d diag(0.5,0.5,0.5);
      diag.normalize();
      Eigen::Vector3d cols(0.1,0.1,0.1);
      int ind = indexfield_(x, y);
      if (ind == -1)
      {
        pixels(x, y) = Col(0);
        continue;
      }
      std::vector<int> inds;
      while (ind != -1)
      {
        inds.push_back(ind);
        ind = trees[ind].attaches_to;
      }
      double scale = 0.7;
      for (int i = (int)inds.size()-1; i>=0; i--)
      {
        srand(1 + inds[i]);
        Eigen::Vector3d hue(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0));
        hue -= diag * hue.dot(diag);
        hue.normalize();
        cols += hue * scale;
        cols += diag * 0.09;
        scale /= 2.0;
      }
      cols[0] = std::max(0.0, std::min(cols[0], 1.0));
      cols[1] = std::max(0.0, std::min(cols[1], 1.0));
      cols[2] = std::max(0.0, std::min(cols[2], 1.0));
      Col col;
      col.a = 255;
      col.r = (uint8_t)(cols[0]*255.0);
      col.g = (uint8_t)(cols[1]*255.0);
      col.b = (uint8_t)(cols[2]*255.0);
      pixels(x, y) = col;
    }
  }

  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawHeightField(const std::string &filename, const Eigen::ArrayXXd &heightfield)
{
  if (!verbose)
    return;
  double max_height = -1e10; 
  double min_height = 1e10;
  for (int i = 0; i<heightfield.rows(); i++)
  {
    for (int j = 0; j<heightfield.cols(); j++)
    {
      if (heightfield(i,j) > -10000.0)
        min_height = std::min(min_height, heightfield(i,j));
      if (heightfield(i,j) < 10000.0)
        max_height = std::max(max_height, heightfield(i,j));
    }
  }

  Field2D<Col> pixels((int)heightfield.rows(), (int)heightfield.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
    for (int y = 0; y < pixels.dims[1]; y++)
      pixels(x, y) = Col((uint8_t)(255.0 * (heightfield(x, y) - min_height)/(max_height - min_height)));
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawGraph(const std::string &filename, const std::vector<Vector4d> &data, double x_min, double x_max, double y_max, double strength_max, double a, double b)
{
  if (!verbose)
    return;
  int res = 256;
  Field2D<Col> pixels(res, res);
  for (auto &c: pixels.data)
    c = Col(0);
  Eigen::Vector3d cols[] = {Eigen::Vector3d(0,0,0), Eigen::Vector3d(1,1,1), Eigen::Vector3d(1,1,0), Eigen::Vector3d(1,0,1), Eigen::Vector3d(0,1,1)};
  for (auto &item: data)
  {
    double x = (double)(res - 1) * (item[0] - x_min) / (x_max - x_min);
    double y = (double)(res - 1) * item[1] / y_max;
    double val = 255.0 * item[2]/strength_max;
    if ((int)item[3] > 1)
      val = 255;
    Eigen::Vector3d c = cols[(int)item[3]];
    if (x >= 0.0 && x<(double)res-1.0 && y >= 0.0 && y<(double)res-1.0)
      pixels((int)x, (int)y) = Col(uint8_t(val*c[0]), uint8_t(val*c[1]), uint8_t(val*c[2]), 255);
  }
  // now draw the line of best fit
  for (int i = 0; i<pixels.dims[0]; i++)
  {
    double x = ((double)i / (double)(res-1)) * (x_max - x_min) + x_min;
    double y = a * (double)x + b;
    int j = (int)((double)(res-1) * y / y_max);
    if (j >= 0 && j < res-1)
      pixels(i, j) += Col(0,127,0,255);
  }

  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawTrees(const std::string &filename, const std::vector<Forest::Result> &results, int width, int height)
{
  double max_height = 0.0;
  double min_height = 1e10;
  for (auto &res: results)
  {
    max_height = std::max(max_height, res.tree_tip[2]);
    min_height = std::min(min_height, res.ground_height);
  }

  // I should probably draw the result
  if (!verbose)
    return;
  Field2D<Col> pixels(width, height);
  for (auto &c: pixels.data)
    c = Col(0); 
  for (auto &result: results)
  {
    Eigen::Vector3d pos = result.tree_tip;
    pos[0] /= voxel_width_;
    pos[1] /= voxel_width_;
 //   double length = pos[2] - result.ground_height;
    double curvature = result.curvature;
    double draw_radius = result.radius; // std::min(0.9 * crown_radius, 50.0); 
    for (int x = (int)(pos[0] - draw_radius); x<= (int)(pos[0]+draw_radius); x++)
    {
      for (int y = (int)(pos[1] - draw_radius); y<= (int)(pos[1]+draw_radius); y++)
      {
        if (x < 0 || x >= width || y<0 || y>=height)
          continue;
        double X = (double)x - pos[0];
        double Y = (double)y - pos[1];
        double mag2 = (double)(X*X + Y*Y);
        if (mag2 <= draw_radius*draw_radius)
        {
          double height = pos[2] + mag2 * curvature;
          double shade = (height - min_height)/(max_height - min_height);
          if (shade > 1.0001)
            std::cout << "weird ass, h: " << height << " pos[2]: " << pos[2] << ", curv: " << curvature << ", mag2: " << mag2 << ", max: " << max_height << std::endl;
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