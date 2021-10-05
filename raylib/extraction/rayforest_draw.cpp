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
// Colour structure
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

// A 2D array of colours
struct ColourField
{
  ColourField() {}
  ColourField(int x, int y){ init(Eigen::Vector2i(x, y)); }
  ColourField(const Eigen::Vector2i &dimensions){ init(dimensions); }
  inline void init(const Eigen::Vector2i &dimensions)
  {
    dims = dimensions;
    data.resize(dims[0] * dims[1]);
  }
  inline Col &operator()(const Eigen::Vector2i &ind){ return data[ind[0] + dims[0]*ind[1]]; } 
  inline const Col &operator()(const Eigen::Vector2i &ind) const { return data[ind[0] + dims[0]*ind[1]]; }
  inline Col &operator()(int x, int y){ return data[x + dims[0]*y]; } 
  inline const Col &operator()(int x, int y) const { return data[x + dims[0]*y]; }
  std::vector<Col> data;
  Eigen::Vector2i dims;
};

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

  ColourField pixels((int)heightfield.rows(), (int)heightfield.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
    for (int y = 0; y < pixels.dims[1]; y++)
      pixels(x, y) = Col((uint8_t)(255.0 * (heightfield(x, y) - min_height)/(max_height - min_height)));
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Occupancy2D::draw(const std::string &filename)
{
  ColourField pixels(dims_[0], dims_[1]);
  for (int x = 0; x < pixels.dims[0]; x++)
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      double shade = pixel(Eigen::Vector3i(x, y, 0)).density();
      pixels(x, y) = Col((uint8_t)(255.0 * shade));
    }
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawFinalSegmentation(const std::string &filename, std::vector<TreeNode> &trees, std::vector<int> &indices)
{
  if (!verbose)
    return;
  ColourField pixels((int)indexfield_.rows(), (int)indexfield_.cols());
  for (int x = 0; x < pixels.dims[0]; x++)
  {
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      int ind = indexfield_(x, y);
      Col col;
      if (ind == -1)
        pixels(x, y) = Col(30);
      else
      {
        while (trees[ind].attaches_to != -1)
        {
          if (std::find(indices.begin(), indices.end(), ind) != indices.end())
            break;
          ind = trees[ind].attaches_to;
        }
        if (std::find(indices.begin(), indices.end(), ind) == indices.end())
        {
          col.a = 255;
          col.r = 255;
          col.g = 0;
          col.b = 255;
          pixels(x, y) = col;
          continue;
        }
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
  ColourField pixels((int)indexfield_.rows(), (int)indexfield_.cols());
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
      double scale = 0.015;
      for (int i = (int)inds.size()-1; i>=0; i--)
      {
        srand(1 + inds[i]);
        Eigen::Vector3d hue(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0));
        hue -= diag * hue.dot(diag);
        hue.normalize();
        cols += hue * scale;
        cols += diag * 0.01;
    //    scale /= 1.25;
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
/*
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
  ColourField pixels(width, height);
  for (auto &c: pixels.data)
    c = Col(0); 
  for (auto &result: results)
  {
    Eigen::Vector3d pos = result.tree_tip;
 //   double length = pos[2] - result.ground_height;
    double curvature = result.curvature;
    double radius_pixels = result.radius / voxel_width_;
    for (int x = (int)(pos[0] - radius_pixels); x<= (int)(pos[0]+radius_pixels); x++)
    {
      for (int y = (int)(pos[1] - radius_pixels); y<= (int)(pos[1]+radius_pixels); y++)
      {
        if (x < 0 || x >= width || y<0 || y>=height)
          continue;
        double X = ((double)x - pos[0]) * voxel_width_;
        double Y = ((double)y - pos[1]) * voxel_width_;
        double mag2 = (double)(X*X + Y*Y);
        if (mag2 <= result.radius*result.radius)
        {
          double height = pos[2] + mag2 * curvature;
          double shade = std::min((height - min_height)/(max_height - min_height), 1.0); // clamp because curvature can conceivably negative sometimes
          Col col(uint8_t(255.0*shade));
          if (pixels(x, y).r < col.r)
            pixels(x, y) = col;
        }
      }
    }
  }    
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}*/

} // namespace ray