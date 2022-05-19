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
#include "raylib/raycloudwriter.h"
#include "raytrees.h"

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
  stbi_flip_vertically_on_write(1);
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
  stbi_flip_vertically_on_write(1);
  stbi_write_png(filename.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void segmentCloud(const std::string &cloud_name_stub, const ColourField &pixels, const Eigen::Vector3d &min_bounds, double voxel_width)
{
  std::string filename = cloud_name_stub + ".ply";
  ray::CloudWriter writer;
  if (!writer.begin(cloud_name_stub + "_segmented.ply"))
    return;

  ray::Cloud chunk;
  auto segment = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    chunk.resize(ends.size());    
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d ind = (ends[i] - min_bounds)/voxel_width;
      chunk.starts[i] = starts[i];
      chunk.ends[i] = ends[i];
      chunk.times[i] = times[i];
      RGBA col;
      col.alpha = colours[i].alpha;
      Col pix = pixels((int)ind[0], (int)ind[1]);
      col.red = pix.r;
      col.green = pix.g;
      col.blue = pix.b;
      chunk.colours[i] = col;
    }
    writer.writeChunk(chunk);
  };

  if (!ray::Cloud::read(filename, segment))
    return;  
  writer.end();
}


void Forest::drawFinalSegmentation(const std::string &cloud_name_stub, std::vector<TreeNode> &trees)
{
  ColourField pixels((int)indexfield_.rows(), (int)indexfield_.cols());
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
          ind = trees[ind].attaches_to;
        if (trees[ind].area <= min_area_)
        {
          pixels(x, y) = Col(0);
        }
        else
        {
          RGBA colour;
          col.a = 255;
          convertIntToColour(ind, colour);
          col.r = colour.red;
          col.g = colour.green;
          col.b = colour.blue;
          pixels(x, y) = col;
        }
      }
    }
  }
  stbi_flip_vertically_on_write(1);

  segmentCloud(cloud_name_stub, pixels, min_bounds_, voxel_width_);

  std::string output_file = cloud_name_stub + "_segmented.png";
  
  stbi_write_png(output_file.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}

void Forest::drawAgglomeration(const std::vector<Cluster> &point_clusters, const std::vector<Eigen::Vector3d> &verts, const std::string &cloud_name_stub)
{
  if (!verbose)
    return;
  ColourField pixels((int)heightfield_.rows(), (int)heightfield_.cols());

  std::vector<double> times;
  Col colour;
  colour.r = colour.g = colour.b = 0;
  colour.a = 255;
  for (int x = 0; x < pixels.dims[0]; x++)
  {
    for (int y = 0; y < pixels.dims[1]; y++)
    {
      pixels(x, y) = colour;
    }
  }  
  for (auto &cluster: point_clusters)
  {
    if (cluster.ids.empty())
      continue;
    colour.r = uint8_t(rand()%255);
    colour.g = uint8_t(rand()%255);
    colour.b = uint8_t(rand()%255);
    for (auto &i: cluster.ids)
    {
      Eigen::Vector3i coord = (verts[i] / voxel_width_).cast<int>();
      for (int x = std::max(0, coord[0] - 1); x <= std::min(coord[0] + 1, pixels.dims[0]-1); x++)
      {
        for (int y = std::max(0, coord[1] - 1); y <= std::min(coord[1] + 1, pixels.dims[1]-1); y++)
        {
          pixels(x, y) = colour;
        }
      }
    }
  }
  stbi_flip_vertically_on_write(1);
  std::string output_file = cloud_name_stub + "_agglom_trees.png";
  stbi_write_png(output_file.c_str(), pixels.dims[0], pixels.dims[1], 4, (void *)&pixels.data[0], 4 * pixels.dims[0]);
}



} // namespace ray