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
  void operator +=(const Col &col)
  {
    r = (uint8_t)std::min((int)r + (int)col.r, 255);
    g = (uint8_t)std::min((int)g + (int)col.g, 255);
    b = (uint8_t)std::min((int)b + (int)col.b, 255);
    a = (uint8_t)std::min((int)a + (int)col.a, 255);
  }
  uint8_t r, g, b, a;
};

static const int res = 256;
static const double max_tree_height = (double)res / 8.0;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;

struct TreeNode
{
  TreeNode() : min_bound(1e10,1e10), max_bound(-1e10,-1e10), attaches_to(-1) 
  {
    curv_mat.setZero();
    curv_vec.setZero();
    sum_square_residual = 0.0;
    sum_square_total = 0.0;
    children[0] = children[1] = -1;
    peak.setZero();
  }
  TreeNode(int x, int y, double height_)
  {
    curv_mat.setZero();
    curv_vec.setZero();
    attaches_to = -1;
    min_bound = max_bound = Eigen::Vector2i(x,y);
    addSample(x,y,height_);
    sum_square_residual = 0.0;
    sum_square_total = 0.0;
    children[0] = children[1] = -1;
    peak = Eigen::Vector3d((double)x, (double)y, height_);
  }
  // for calculating paraboloid of best fit:
  Matrix4d curv_mat;
  Vector4d curv_vec;
  Vector4d abcd; // the solved paraboloid
  double sum_square_residual;
  double sum_square_total;
  Eigen::Vector2i min_bound, max_bound;
  Eigen::Vector3d peak;
  int attaches_to;
  int children[2];

  Eigen::Vector2d centroid() const { return Eigen::Vector2d(curv_mat(1,3) / area(), curv_mat(2,3) / area()); }
  inline double area() const { return curv_mat(3,3); }
  inline double avgHeight() const { return curv_vec[3] / area(); }
  inline double height() const { return abcd[3] - (abcd[1]*abcd[1] + abcd[2]*abcd[2])/(4*abcd[0]); }
  inline Eigen::Vector3d tip() const { return Eigen::Vector3d(-abcd[1]/(2*abcd[0]), -abcd[2]/(2*abcd[0]), height()); }
  inline double heightAt(double x, double y) const { return abcd[0]*(x*x + y*y) + abcd[1]*x + abcd[2]*y + abcd[3]; }
  inline double crownRadius() const { return 1.0 / -abcd[0]; }
  inline bool validParaboloid() const 
  {
    const double minimum_crown_radius = 0.5;
    const double maximum_crown_radius = 20.0; // in metres
    double r = crownRadius();
    if (r<minimum_crown_radius || r > maximum_crown_radius)
      return false;
    Eigen::Vector3d top = tip();
    for (int i = 0; i<2; i++)
      if (top[i] < (double)min_bound[i] || top[i] > (double)max_bound[i])
        return false;
    return true;
  }
  inline void addSample(double x, double y, double z)
  {
    Vector4d vec(x*x + y*y, x, y, 1.0);
    curv_mat += vec * vec.transpose();
    curv_vec += z*vec;
  }
  void updateBound(const Eigen::Vector2i &bmin, const Eigen::Vector2i &bmax)
  {
    for (int i = 0; i<2; i++)
    {
      min_bound[i] = std::min(min_bound[i], bmin[i]);
      max_bound[i] = std::max(max_bound[i], bmax[i]);
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
        srand(1 + ind);
        Col col;
        col.a = 255;
        col.r = (uint8_t)(rand()%256);
        col.g = (uint8_t)(rand()%256);
        col.b = (uint8_t)(rand()%256);
        pixels[x + res * y] = col;
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
    if (X>=0 && X<res && Y >= 0.0 && Y<res)
      pixels[X + res*Y] = Col((uint8_t)(255.0 * ray::clamped(z/max_tree_height, 0.0, 1.0)));
  }

  stbi_write_png("segmenting.png", res, res, 4, (void *)&pixels[0], 4 * res);
}

void searchTrees(const std::vector<TreeNode> &trees, int ind, double error, double length_per_radius, double ground_height, std::vector<int> &indices)
{
  if (trees[ind].children[0] == -1)
  {
    if (trees[ind].validParaboloid())
      indices.push_back(ind);
    return;
  }
  int ind0 = trees[ind].children[0];
  double base0 = trees[ind0].height() - length_per_radius * trees[ind0].crownRadius();
  double error0 = abs(base0 - ground_height);
  int ind1 = trees[ind].children[1];
  double base1 = trees[ind1].height() - length_per_radius * trees[ind1].crownRadius();
  double error1 = abs(base1 - ground_height);
      
  if (error < std::min(error0, error1) && trees[ind].validParaboloid()) // we've found the closest, so end loop
  {
    indices.push_back(ind);
    return;
  }
  searchTrees(trees, ind0, error0, length_per_radius, ground_height, indices);
  searchTrees(trees, ind1, error1, length_per_radius, ground_height, indices);
}

// Decimates the ray cloud, spatially or in time
int main(int /*argc*/, char */*argv*/[])
{
  const bool verbose = true;
  for (int seed = 0; seed<=300; seed+=100)
  {
  double heightfield[res][res];
  int indexfield[res][res];
  for (int i = 0; i<res; i++)
    for (int j = 0; j<res; j++)
      indexfield[i][j] = -1;
  memset(heightfield, 0, res*res*sizeof(double));
  #if 1 // create the height field
  int num = 500;
  const double radius_to_height = 0.4;
  if (seed)
    srand(seed);
//  srand(200);
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
  if (1)
  {
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
      }
    }
  }
  double max_height = 0.0;
  for (int i = 0; i<res; i++)
    for (int j = 0; j<res; j++)
      max_height = std::max(max_height, heightfield[i][j]);
  // now render it 
  if (verbose)
  {
    std::vector<Col> pixels(res * res);
    for (int x = 0; x < res; x++)
      for (int y = 0; y < res; y++)
        pixels[x + res * y] = Col((uint8_t)(255.0 * heightfield[x][y]/max_height));
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
    int x, y, index; // if index == -2 then we are merging
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
  std::set<int> heads;
  // 1. find highest points
  for (int x = 0; x < res; x++)
  {
    for (int y = 0; y < res; y++)
    {
      // Moore neighbourhood
      double height = heightfield[x][y];
      double max_h = 0.0;
      for (int i = std::max(0, x-1); i<= std::min(x+1, res-1); i++)
        for (int j = std::max(0, y-1); j<= std::min(y+1, res-1); j++)
          if (!(i==x && j==y))
            max_h = std::max(max_h, heightfield[i][j]);
      if (height > max_h)
      {
        Point p;
        p.x = x; p.y = y; p.height = height;
        p.index = (int)basins.size();
        basins.insert(p);
        heads.insert((int)trees.size());
        trees.push_back(TreeNode(x,y,height));
      }
    }
  }
  std::cout << "initial number of peaks: " << trees.size() << std::endl;
  // now iterate until basins is empty
  int cnt = 0;
  const int max_tree_length = 22;
  while (!basins.empty())
  {
    Point p = *basins.begin();
    int x = p.x;
    int y = p.y;
    basins.erase(p); // removes it from basins. p still exists

    if (p.index == -2) // a merge request
    {
      int p_head = x;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
      int q_head = y;
      while (trees[q_head].attaches_to != -1)
        q_head = trees[q_head].attaches_to;
      if (p_head != q_head)
      {
        TreeNode &p_tree = trees[p_head];
        TreeNode &q_tree = trees[q_head];
        Eigen::Vector2i mx = ray::maxVector2(p_tree.max_bound, q_tree.max_bound);
        Eigen::Vector2i mn = ray::minVector2(p_tree.min_bound, q_tree.min_bound);
        mx -= mn;
        if (std::max(mx[0], mx[1]) <= max_tree_length)
        {
          int new_index = (int)trees.size();
          TreeNode node;
          node.curv_mat = p_tree.curv_mat + q_tree.curv_mat;
          node.curv_vec = p_tree.curv_vec + q_tree.curv_vec;
          node.min_bound = p_tree.min_bound;
          node.max_bound = p_tree.max_bound;
          node.updateBound(q_tree.min_bound, q_tree.max_bound);
          node.children[0] = p_head;  
          node.children[1] = q_head;
          node.abcd = node.curv_mat.ldlt().solve(node.curv_vec);

          if (node.validParaboloid())
          {
            heads.erase(p_head);
            heads.erase(q_head);
            heads.insert(new_index);
            p_tree.attaches_to = new_index;
            q_tree.attaches_to = new_index;
            trees.push_back(node); // danger, this can invalidate the p_tree reference
          }
        }
      }
      continue;
    }    
    indexfield[x][y] = p.index;

    int xs[4] = {x-1, x, x, x+1};
    int ys[4] = {y, y+1, y-1, y};
    for (int i = 0; i<4; i++)
    {
      if (xs[i] < 0 || xs[i] >= res)
        continue;
      if (ys[i] < 0 || ys[i] >= res)
        continue;
      int p_head = p.index;
      while (trees[p_head].attaches_to != -1)
        p_head = trees[p_head].attaches_to;
        
      int xx = xs[i];
      int yy = ys[i];
      int &ind = indexfield[xx][yy];

      int q_head = ind;
      if (q_head != -1)
      {
        while (trees[q_head].attaches_to != -1)
          q_head = trees[q_head].attaches_to;
      }

      if (ind != -1 && p_head != q_head)
      {
        TreeNode &p_tree = trees[p_head];
        TreeNode &q_tree = trees[q_head];
        cnt++;
//        if (verbose && !(cnt%50)) // I need a way to visualise the hierarchy here!
//          drawSegmentation(indexfield, trees);
        Eigen::Vector2i mx = ray::maxVector2(p_tree.max_bound, q_tree.max_bound);
        Eigen::Vector2i mn = ray::minVector2(p_tree.min_bound, q_tree.min_bound);
        mx -= mn;
        bool merge = std::max(mx[0], mx[1]) <= max_tree_length;
        if (merge)
        {
          const double flood_merge_scale = 1.5; // 1 merges immediately, infinity never merges
          // add a merge task:
          Eigen::Vector2d mid(xx, yy);
          double p_sqr = (Eigen::Vector2d(p_tree.peak[0], p_tree.peak[1]) - mid).squaredNorm();
          double q_sqr = (Eigen::Vector2d(q_tree.peak[0], q_tree.peak[1]) - mid).squaredNorm();
          double blend = p_sqr / (p_sqr + q_sqr);
          double flood_base = p_tree.peak[2]*(1.0-blend) + q_tree.peak[2]*blend;
          double low_flood_height = flood_base - p.height;

          Point q;
          q.x = p_head; q.y = q_head; 
          q.index = -2;
          q.height = flood_base - low_flood_height * flood_merge_scale;
          basins.insert(q);
        }
      }
      if (ind == -1 && heightfield[xx][yy] > 0.0)
      {
  //      Eigen::Vector2i dif = (trees[p_head].max_bound - trees[p_head].min_bound);
  //      double len = 1.0 + (dif[0] + dif[1])/2.0;
  //      double height_delta = len * 0.5;
        Point q;
        q.x = xx; q.y = yy; q.index = p.index;
        q.height = heightfield[xx][yy];
   //     if ((p.height - q.height) < 3.0 + len * 0.5)
   //     {
        ind = p.index;
        basins.insert(q);
        trees[p_head].addSample(xx, yy, q.height);
        trees[p_head].updateBound(Eigen::Vector2i(xx, yy), Eigen::Vector2i(xx, yy));
  //      }
      }
    }
  }
  if (verbose) // I need a way to visualise the hierarchy here!
  {
    drawSegmentation(indexfield, trees);
    std::vector<Col> pixels(res * res);
    for (int x = 0; x < res; x++)
    {
      for (int y = 0; y < res; y++)
      {
        srand(1 + indexfield[x][y]);
        Col col;
        col.a = 255;
        col.r = uint8_t(rand()%256);
        col.g = uint8_t(rand()%256);
        col.b = uint8_t(rand()%256);
        pixels[x + res * y] = indexfield[x][y] == -1 ? Col(0) : col;
      }
    }
    stbi_write_png("segmented.png", res, res, 4, (void *)&pixels[0], 4 * res);
  }

//  std::vector<Node> nodes(attachto.size());
  std::cout << "number of raw candidates: " << trees.size() << " number largest size: " << heads.size() << std::endl;

  // calculate features of leaf nodes
  double max_area = 0.0;
  for (auto &tree: trees)
    max_area = std::max(max_area, tree.area());
  // now plot curvatures against heights, to find a correlation..

  // first I'll have to solve each tree parabola and store it in a vector:
  for (auto &tree: trees)
    tree.abcd = tree.curv_mat.ldlt().solve(tree.curv_vec);

  // for each pixel, we have to update accuracy data on each tree.
  for (int x = 0; x < res; x++)
  {
    for (int y = 0; y < res; y++)
    {
      int ind = indexfield[x][y];
      if (ind != -1)
      {
        TreeNode &tree = trees[ind];
        tree.sum_square_residual += ray::sqr(heightfield[x][y] - tree.heightAt(x,y));
        tree.sum_square_total += ray::sqr(heightfield[x][y] - tree.avgHeight());
        while (trees[ind].attaches_to != -1)
        {
          ind = trees[ind].attaches_to;
          TreeNode &tree = trees[ind];
          tree.sum_square_residual += ray::sqr(heightfield[x][y] - tree.heightAt(x,y));
          tree.sum_square_total += ray::sqr(heightfield[x][y] - tree.avgHeight());
        }
      }
    }
  }

  std::vector<Col> pixels(res * res);
  std::vector<Col> pixels2(res * res);
  for (auto &c: pixels)
    c = Col(20);
  for (auto &c: pixels2)
    c = Col(0);
  // so radius = height * radius_to_height
  // parabola raises height/2 in radius, so height/2 = curvature*(radius^2)
  // so predicted height = rad/(2*radius_to_height^2)
  std::vector<Eigen::Vector3d> data;
  // draw scatter plot
  for (auto &tree: trees)
  {
    if (!tree.validParaboloid())
      continue;
    double x = (double)(res - 1) * tree.height() / max_tree_height;
    double predicted_length = (tree.crownRadius() / (2.0 * ray::sqr(radius_to_height)));
    double y = (double)(res - 1) * predicted_length / (2.0 * max_tree_height);
    double strength = 255.0 * std::sqrt(tree.area() / 300.0);
    double strength_sqr = 255.0 * (tree.area() / 300.0);
    double ratio = tree.sum_square_residual/(1e-10 + tree.sum_square_total);
    double R2 = 1.0 - std::sqrt(ray::clamped(ratio, 0.0, 1.0)); 
    strength *= R2; // TODO: whether this helps is dubious, and it is weird that it gives values outside 0-1... more testing needed.
    if (x==x && x>=0.0 && x<(double)res)
    {
      if (y==y && y>=0.0 && y<(double)res)
      { 
        data.push_back(Eigen::Vector3d(tree.height(), tree.crownRadius(), strength));
        pixels[(int)x + res*(int)y] += Col((uint8_t)std::min(strength, 255.0));
      }

      double y2 = (double)(res - 1) * sqrt(tree.area() / max_area) * 0.5;
      if (y2 >= 0.0 && y2 < (double)res)
        pixels2[(int)x + res*(int)y2] += Col((uint8_t)std::min(strength_sqr, 255.0));
    }
  }
  // now analyse the data to get line of best fit:
  // line = radius = a * treetop_pos + b
  double a = 0.0;
  double b = 0.0; // Nope! this is the wrong way around, so this is not ground_height!
  double min_height_resolving = 0.5; // unlikely to get better resolution than this
  for (int it = 0; it<13; it++) // from mean to median to mode (ish)
  {
    double power = (double)it / 9.0; // just 1 is median... its very similar, but not obviously better
    Eigen::Vector3d mean(0,0,0);
    double total_weight = 0.0;
    for (auto &point: data)
    {
      double error = abs(point[1] - (a * point[0] + b));
      double weight = point[2] / (min_height_resolving + pow(error, power));
      mean += point * weight;
      total_weight += weight;
    }
    mean /= total_weight;
    double total_xy = 0.0;
    double total_xx = 0.0;
    for (auto &point: data)
    {
      double error = abs(point[1] - (a * point[0] + b));
      double weight = point[2] / (min_height_resolving + pow(error, power));
      Eigen::Vector3d p = point - mean;
      total_xy += p[0]*p[1] * weight;
      total_xx += p[0]*p[0] * weight;
    }
    a = total_xy / std::max(1e-10, total_xx);
    b = mean[1] - a*mean[0];
  }
  std::cout << "a: " << a << ", b: " << b << std::endl;

  // draw it:
 /* for (int x = 0; x<res; x++)
  {
    double y = gradient * (double)x + ground_height;
    int Y = (int)y;
    if (Y >=0 && Y < res)
      pixels[x + res*Y] += Col(64);
  }*/

  stbi_write_png("graph_curv.png", res, res, 4, (void *)&pixels[0], 4 * res);
  stbi_write_png("graph_area.png", res, res, 4, (void *)&pixels2[0], 4 * res);
  // radius = a*treetop_pos + b, so treetop_pos = (radius - b)/a
  double ground_height = -b/a;
  double length_per_radius = 1.0/a;

  /***************************************************************************************/
  // The last step is to take the ground surface (currently flat: 'b') and scale estimation 'a'
  // and use this to estimate the set of tree heights
  // Note: we may still choose an extrapolation technique later... this is worth pursuing at some point
  struct Result
  {
    std::vector<Eigen::Vector3d> tree_tips;
    double ground_height;
    double treelength_per_crownradius;
  };
  Result result;
  max_height = ground_height;
  std::vector<int> indices;
  for (auto &head: heads)
  {
    int ind = head;
    double base = trees[ind].height() - length_per_radius * trees[ind].crownRadius();
    double error = abs(base - ground_height);
    searchTrees(trees, head, error, length_per_radius, ground_height, indices);
  }
  for (auto &ind: indices)
  {
    result.tree_tips.push_back(trees[ind].tip());
    max_height = std::max(max_height, trees[ind].height());
  }
  max_height = max_tree_height * 1.2; // TODO: fix!, problem is one of the tree heights is really big, 2 lines up
  result.ground_height = ground_height;
  result.treelength_per_crownradius = length_per_radius;
  // I should probably draw the result
  if (verbose)
  {
    std::vector<Col> pixels(res * res);
    for (auto &c: pixels)
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
            if (pixels[x + res*y].r < col.r)
              pixels[x + res*y] = col;
          }
        }
      }
    }    
    stbi_write_png("result_trees.png", res, res, 4, (void *)&pixels[0], 4 * res);
  }
  } 
  return 0;
}
