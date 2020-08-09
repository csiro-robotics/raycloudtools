// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayutils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "raylib/imagewrite.h"

#define _TEST // run multiple times

void usage(bool error=false)
{
  std::cout << "Extract feature into a text file structure" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayextract trees method 1    - for low coverage scans (e.g. flyovers), where the trees are approximately inferred" << std::endl;
//  cout << " --extrapolate  - estimates tree distribution and adds trees where there is no evidence to the contrary" << endl;
  exit(error);
}

struct Col
{
  Col(){}
  Col(uint8_t shade) : r(shade), g(shade), b(shade), a(255) {}
  uint8_t r, g, b, a;
};

static const int res = 256;

struct TreeNode
{
  TreeNode() : centroid(0,0), height(0), area(1), curvature(0), min_bound(1e10,1e10), max_bound(-1e10,-1e10), attaches_to(-1) {}
  Eigen::Vector2d centroid;
  double height;
  double area;
  double curvature;
  Eigen::Vector2i min_bound, max_bound;
  int attaches_to;
  int length() const
  {
    Eigen::Vector2i dif = max_bound - min_bound;
    if (dif[0] < 0)
      dif[0] += res;
    if (dif[1] < 0)
      dif[1] += res;
    return std::max(dif[0], dif[1]);
  }
  void updateBound(const Eigen::Vector2i &bmin, const Eigen::Vector2i &bmax)
  {
    for (int i = 0; i<2; i++)
    {
      int diff = (res + (bmax[i] - max_bound[i]))%res;
      if (diff < res/2)
        max_bound[i] = bmax[i];

      int diff2 = (res + (min_bound[i] - bmin[i]))%res;
      if (diff2 < res/2)
        min_bound[i] = bmin[i];
    }
  }
};

void drawSegmentation(int indexfield[][res], const std::vector<TreeNode> &trees)
{
  std::vector<Col> pixels(res * res);
  for (int x = 0; x < res; x++)
  {
    for (int y = 0; y < res; y++)
    {
      int ind = indexfield[x][y];
      if (ind == -1)
        pixels[x + res * y] = Col(0);
      else
      {
        while (trees[ind].attaches_to != -1)
          ind = trees[ind].attaches_to;
        if (0)
        {
          Eigen::Vector2i d = trees[ind].max_bound - trees[ind].min_bound;
          if (d[0] >= 0)
          {
            if (x <trees[ind].min_bound[0] || x > trees[ind].max_bound[0])
              std::cout <<"bad bounds" << std::endl;
          }
          else
          {
            if (x < trees[ind].min_bound[0] && x > trees[ind].max_bound[0])
              std::cout <<"bad bounds" << std::endl;
          }
          if (d[1] >= 0)
          {
            if (y <trees[ind].min_bound[1] || y > trees[ind].max_bound[1])
              std::cout <<"bad bounds" << std::endl;
          }
          else
          {
            if (y < trees[ind].min_bound[1] && y > trees[ind].max_bound[1])
              std::cout <<"bad bounds" << std::endl;
          }
        }

        srand(1 + ind);
        Col col;
        col.a = 255;
        col.r = rand()%256;
        col.g = rand()%256;
        col.b = rand()%256;
        pixels[x + res * y] = col;
      }
    }
  }
  /*for (int i = 0; i<(int)trees.size(); i++)
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
  stbi_write_png("segmenting.png", res, res, 4, (void *)&pixels[0], 4 * res);
}

// Decimates the ray cloud, spatially or in time
int main(int argc, char *argv[])
{
  const bool verbose = true;
  double heightfield[res][res];
  int indexfield[res][res];
  for (int i = 0; i<res; i++)
    for (int j = 0; j<res; j++)
      indexfield[i][j] = -1;
  memset(heightfield, 0, res*res*sizeof(double));
  #if 1 // create the height field
  int num = 500;
  double max_tree_height = (double)res / 8.0;
  const double radius_to_height = 0.4;
  std::vector<Eigen::Vector3d> ps(num);
  for (int i = 0; i<num; i++)
  {
    double height = max_tree_height * std::pow(2.0, ray::random(-2.0, 0.0)); // i.e. 0.25 to 1 of max_height
    ps[i] = Eigen::Vector3d(ray::random(0.0, (double)res-1.0), ray::random(0.0, (double)res - 1.0), height);
  }
  for (int it = 0; it < 10; it++)
  {
    for (auto &p: ps)
    {
      Eigen::Vector3d shift(0,0,0);
      for (auto &other: ps)
      {
        double radius = radius_to_height * (p[2] + other[2]);
        Eigen::Vector3d dif = p - other;
        dif[2] = 0.0;
        if (dif[0] > res/2.0)
          dif[0] -= res;
        if (dif[1] > res/2.0)
          dif[1] -= res;
        double len_sqr = dif.squaredNorm();
        if (len_sqr > 0.0 && len_sqr < radius*radius)
        {
          double len = std::sqrt(len_sqr);
          shift += (dif/len) * (radius-len);
        }
      }
      p += 0.5*shift;
      p[0] = fmod(p[0] + (double)res, (double)res);
      p[1] = fmod(p[1] + (double)res, (double)res);
    }
  }
  // now make height field
  double max_height = 0.0;
  for (auto &p: ps)
  {
    double radius = radius_to_height * p[2];
    for (int x = (int)(p[0] - radius); x<= (int)(p[0]+radius); x++)
    {
      for (int y = (int)(p[1] - radius); y<= (int)(p[1]+radius); y++)
      {
        double X = (double)x - p[0];
        double Y = (double)y - p[1];
        double mag2 = (double)(X*X + Y*Y);
        if (mag2 <= radius*radius)
        {
          double height = p[2] - (p[2]/2.0)*mag2/(radius*radius);
          int xx = (x + res)%res;
          int yy = (y + res)%res;
          height += ray::random(-1.0, 1.0);
          heightfield[xx][yy] = std::max(heightfield[xx][yy], height);
        }
      }
    }
  }
  // add a noisy function:
  const int wid = 80;
  double hs[wid][wid];
  for (int i = 0; i<wid; i++)
  {
    for (int j = 0; j<wid; j++)
      hs[i][j] = ray::random(-5.0, 5.0);
  }
  for (int i = 0; i<res; i++)
  {
    double x = (double)wid * (double)i/(double)res;
    int X = (int)x;
    double blendX = x-(double)X;
    for (int j = 0; j<res; j++)
    {
      double y = (double)wid * (double)j/(double)res;
      int Y = (int)y;
      double blendY = y-(double)Y;
      heightfield[i][j] += hs[X][Y]*(1.0-blendX)*(1.0-blendY) + hs[X][Y+1]*(1.0-blendX)*blendY + 
                           hs[X+1][Y]*blendX*(1.0-blendY) + hs[X+1][Y+1]*blendX*blendY; 
      heightfield[i][j] = std::max(heightfield[i][j], 0.0);
      max_height = std::max(max_height, heightfield[i][j]);
    }
  }
  // now render it 
  if (verbose)
  {
    std::vector<Col> pixels(res * res);
    for (int x = 0; x < res; x++)
      for (int y = 0; y < res; y++)
        pixels[x + res * y] = Col(255.0 * heightfield[x][y]/max_height);
    stbi_write_png("testheight.png", res, res, 4, (void *)&pixels[0], 4 * res);
  }
  #endif

  // Now, we want a watershed algorithm to pick the individual shapes:
  // Algorithm: 
  // 1. find all the highest points, give them unique labels, sort them from highest to lowest
  // 2. for each highest point, make it black, then give all the non-black neighbours the same label and 
  //    add them to the sorted list

  struct Point 
  { 
    int x, y, index; 
    double height;
  };
  struct PointCmp 
  {
    bool operator()(const Point& lhs, const Point& rhs) const 
    { 
      return lhs.height > rhs.height; 
    }
  };  
  std::set<Point, PointCmp> basins;
  std::vector<TreeNode> trees;
  // 1. find highest points
  for (int x = 0; x < res; x++)
  {
    for (int y = 0; y < res; y++)
    {
      // Moore neighbourhood
      int xs[8] = {x-1, x, x+1, x-1, x+1, x-1, x, x+1};
      int ys[8] = {y-1, y-1, y-1, y, y, y+1,y+1,y+1};
      double height = heightfield[x][y];
      double max_h = 0.0;
      for (int i = 0; i<8; i++)
        max_h = std::max(max_h, heightfield[(xs[i]+res)%res][(ys[i]+res)%res]);
      if (height > max_h)
      {
        Point p;
        p.x = x; p.y = y; p.height = height;
        p.index = basins.size();
        basins.insert(p);
        TreeNode node;
        node.centroid = Eigen::Vector2d(x,y);
        node.height = height;
        node.min_bound = node.max_bound = Eigen::Vector2i(x,y);
        trees.push_back(node);
      }
    }
  }
  std::cout << "initial number of peaks: " << trees.size() << std::endl;
  // now iterate until basins is empty
  int cnt = 0;
  int max_tree_length = 20;
  while (!basins.empty())
  {
    Point p = *basins.begin();
    int x = p.x;
    int y = p.y;
    basins.erase(p); // how do I erase the first member, it shouldnt need to search for this
    indexfield[x][y] = p.index;

    int xs[4] = {x-1, x, x, x+1};
    int ys[4] = {y, y+1, y-1, y};
    for (int i = 0; i<4; i++)
    {
      int p_head = p.index;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
      TreeNode &p_tree = trees[p_head];
        
      int xx = (xs[i] + res)%res;
      int yy = (ys[i] + res)%res;
      int &ind = indexfield[xx][yy];

      int q_head = ind;
      if (q_head != -1)
      {
        while (trees[q_head].attaches_to != -1)
          q_head = trees[q_head].attaches_to;
      }

      if (ind != -1 && p_head != q_head)
      {
        TreeNode &q_tree = trees[q_head];
        cnt++;
        if (verbose && !(cnt%50)) // I need a way to visualise the hierarchy here!
          drawSegmentation(indexfield, trees);
        if (std::min(p_tree.area, q_tree.area) > std::max(p_tree.area, q_tree.area) * 0.1)
        {
          bool merge = false;
          if (p_tree.area > q_tree.area)
            merge = p_tree.length() <= max_tree_length;
          else
            merge = q_tree.length() <= max_tree_length;
          if (merge)
          {
            int new_index = trees.size();
            TreeNode node;
            node.area = p_tree.area + q_tree.area;
            node.min_bound = p_tree.min_bound;
            node.max_bound = p_tree.max_bound;
            node.updateBound(q_tree.min_bound, q_tree.max_bound);
            trees.push_back(node);
            p_tree.attaches_to = new_index;
            q_tree.attaches_to = new_index;
          }
        }
        else
        {
          if (p_tree.area > q_tree.area)
          {
            q_tree.attaches_to = p_head;
            p_tree.updateBound(q_tree.min_bound, q_tree.max_bound);
          }
          else
          {
            p_tree.attaches_to = q_head;
            q_tree.updateBound(p_tree.min_bound, p_tree.max_bound);
          }
        }
      }
      if (ind == -1 && heightfield[xx][yy] > 0.0)
      {
        Point q;
        q.x = xx; q.y = yy; q.index = p.index;
        q.height = heightfield[xx][yy];
        ind = p.index;
        basins.insert(q);
        p_tree.area++;
        p_tree.updateBound(Eigen::Vector2i(xx, yy), Eigen::Vector2i(xx, yy));
      }
    }
  }
  if (verbose) // I need a way to visualise the hierarchy here!
  {
    std::vector<Col> pixels(res * res);
    for (int x = 0; x < res; x++)
    {
      for (int y = 0; y < res; y++)
      {
        srand(1 + indexfield[x][y]);
        Col col;
        col.a = 255;
        col.r = rand()%256;
        col.g = rand()%256;
        col.b = rand()%256;
        pixels[x + res * y] = indexfield[x][y] == -1 ? Col(0) : col;
      }
    }
    stbi_write_png("segmented.png", res, res, 4, (void *)&pixels[0], 4 * res);
  }

//  std::vector<Node> nodes(attachto.size());
  std::cout << "number of raw candidates: " << trees.size() << std::endl;
#if 0
  // calculate features of leaf nodes
  double min_curvature = 1e10;
  double max_area = 0.0;
  for (int i = 0; i<(int)trees.size(); i++)
  {
    double curv = 0.0;
    double area = 0.0;
    double sum_dist2 = 0.0;
    for (int x = 0; x < res; x++) // can speed this bit up
    {
      for (int y = 0; y < res; y++)
      {
        if (indexfield[x][y] == i)
        {
          double difX = abs(x - trees[i][0]);
          double difY = abs(y - trees[i][1]);
          if (difX > (double)(res/2))
            difX -= res;
          if (difY > (double)(res/2))
            difY -= res;
          int dist2 = difX*difX + difY*difY;
          curv += (trees[i][2] - heightfield[x][y]);
          sum_dist2 += dist2;
          area++;
        }
      }
    }
    nodes[i].curvature = curv/sum_dist2;
    if (nodes[i].curvature != 0.0)
      min_curvature = std::min(min_curvature, nodes[i].curvature);
    nodes[i].area = area;
    nodes[i].centroid = Eigen::Vector2d(trees[i][0], trees[i][1]);
    nodes[i].height = trees[i][2];
    max_area = std::max(max_area, area);
  }
  // now fill in features from root to tip:
  for (size_t i = trees.size(); i<parents.size(); i++)
  {
    Node n0 = nodes[parents[i][0]];
    Node n1 = nodes[parents[i][1]];
    Node &n = nodes[i];
    n.area = n0.area + n1.area;
    n.centroid = (n0.centroid*n0.area + n1.centroid*n1.area)/(n0.area + n1.area); // Could do better
    n.height = std::max(n0.height, n1.height);
    n.curvature = n0.curvature*n0.area + n1.curvature*n1.area; // this is wrong!
  }
  // now plot curvatures against heights, to find a correlation..
  if (verbose)
  {
    std::vector<Col> pixels(res * res);
    for (auto &c: pixels)
      c = Col(0);
    for (auto &node: nodes)
    {
      double x = (double)(res - 1) * node.height / max_tree_height;
      double y = (double)(res - 1) * min_curvature / node.curvature;
      if (node.curvature != 0.0)
        pixels[(int)x + res*(int)y] = Col(255);
      double y2 = (double)(res - 1) * sqrt(node.area / max_area);
      Col col2(128);
      col2.r = 255;
      pixels[(int)x + res*(int)y2] = col2;
    }
    stbi_write_png("graph.png", res, res, 4, (void *)&pixels[0], 4 * res);
  } 
  #endif
  return 0;
}
