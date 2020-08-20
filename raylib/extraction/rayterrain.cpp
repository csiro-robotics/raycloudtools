// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterrain.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayconvexhull.h"

namespace ray
{
struct Node
{
  Node()
  { 
    pos.setZero();
    found = 0;
    for (int i = 0; i<2; i++)
      for (int j = 0; j<2; j++)
        for (int k = 0; k<2; k++)
          dir_ids[i][j][k] = -1;
    min_bounds = Eigen::Vector3d(-1e10,-1e10,-1e10);
    max_bounds = Eigen::Vector3d(1e10,1e10,1e10);
  }
  Vector4d pos;
  Eigen::Vector3d min_bounds, max_bounds;
  int found;
  int dir_ids[2][2][2];

  bool somethingLarger(std::vector<Node> &nodes, const Vector4d &corner)
  {
    Vector4d dif = corner - pos;
    if (dif == Vector4d(0,0,0,0))
      return dir_ids[1][1][1] != -1;
    int i = (int)(dif[0] > 0.0);
    int j = (int)(dif[1] > 0.0);
    int k = (int)(dif[2] > 0.0);
    if (i==1 && j==1 && k==1) // corner is larger, so deactivate current node
      found = 1; 
    else if (i==0 && j==0 && k==0) // corner is smaller, so this node is indeed larger
      return true;
    if (corner[0] >= max_bounds[0] || corner[1] >= max_bounds[1] || corner[2] >= max_bounds[2])
      return false; // corner is out of bounds, so we have not found something larger
    for (int I = i; I<2; I++)
      for (int J = j; J<2; J++)
        for (int K = k; K<2; K++)
          if (dir_ids[I][J][K] != -1)
            if (nodes[dir_ids[I][J][K]].somethingLarger(nodes, corner))
              return true;
    return false;
  }
  bool somethingSmaller(std::vector<Node> &nodes, const Vector4d &corner)
  {
    Vector4d dif = corner - pos;
    if (dif == Vector4d(0,0,0,0))
      return dir_ids[0][0][0] != -1;
    int i = (int)(dif[0] > 0.0);
    int j = (int)(dif[1] > 0.0);
    int k = (int)(dif[2] > 0.0);
    if (i==0 && j==0 && k==0) // corner is larger, so deactivate current node
      found = 1; 
    else if (i==1 && j==1 && k==1) // corner is smaller, so this node is indeed larger
      return true;
    if (corner[0] <= min_bounds[0] || corner[1] <= min_bounds[1] || corner[2] <= min_bounds[2])
      return false; // corner is out of bounds, so we have not found something larger
    for (int I = 0; I<=i; I++)
      for (int J = 0; J<=j; J++)
        for (int K = 0; K<=k; K++)
          if (dir_ids[I][J][K] != -1)
            if (nodes[dir_ids[I][J][K]].somethingSmaller(nodes, corner))
              return true;
    return false;
  }
};

void constructOctalSpacePartition(std::vector<Node> &nodes, std::vector<Vector4d> points)
{
  nodes.resize(points.size());
  int i = 0;
  while (points.size() > 0)
  {
    int ind = rand()%(int)points.size();
    nodes[i++].pos = points[ind];
    points[ind] = points.back(); points.pop_back();
  }
  // n log n on average, for each node, trace the tree to add correct direction ids
  for (size_t n = 1; n<nodes.size(); n++)
  {
    Eigen::Vector4d &pos = nodes[n].pos;
    int head = 0;
    for (;;)
    {
      Vector4d dif = pos - nodes[head].pos;
      int i = (int)(dif[0] > 0.0);
      int j = (int)(dif[1] > 0.0);
      int k = (int)(dif[2] > 0.0);
      int new_head = nodes[head].dir_ids[i][j][k];
      if (new_head == -1)
      {
        nodes[head].dir_ids[i][j][k] = (int)n;
        nodes[n].min_bounds = nodes[head].min_bounds;
        nodes[n].max_bounds = nodes[head].max_bounds;
        int inds[3] = {i,j,k};
        for (int l = 0; l<3; l++)
        {
          double &bound = inds[l]==1 ? nodes[n].min_bounds[l] : nodes[n].max_bounds[l];
          bound = nodes[head].pos[l];
        }
        break;
      }
      head = new_head;
    }
  }
}

void Terrain::cropPointsAbove(std::vector<Vector4d> &points, std::vector<Vector4d> &neg, double width)
{
  std::vector<Vector4d> front = neg;
  double offset = width / std::sqrt(3.0);
  for (auto &v: front)
  {
    v[0] += offset; v[1] += offset; v[2] += offset;
  }
  std::vector<Node> nodes;
  constructOctalSpacePartition(nodes, front);
  Node &root = nodes[0];
  std::vector<int> new_index(points.size());
  std::vector<int> indices;
  indices.reserve(points.size());
  for (size_t i = 0; i<points.size(); i++)
  {
    if (root.somethingLarger(nodes, points[i]))
    {
      new_index[i] = (int)indices.size();
      indices.push_back((int)i);
    }
    else
      new_index[i] = -1;
  }
  for (auto &v: neg)
  {
    int ind = new_index[(int)(v[3])];
    if (ind == -1)
      std::cout << "error, old points should always be visible from higher ones" << std::endl;
    else
      v[3] = (double)ind + 0.5;
  }
  std::vector<Vector4d> new_points;
  new_points.reserve(indices.size());
  for (auto &i: indices)
  {
    new_points.push_back(points[i]);
    new_points.back()[3] = (double)(new_points.size()-1) + 0.5;
  }
  points = new_points;
}


void Terrain::getParetoFront(const std::vector<Vector4d> &points, std::vector<Vector4d> &front, bool positive)
{
  std::vector<Node> nodes;
  constructOctalSpacePartition(nodes, points);
  Node &root = nodes[0];
  for (size_t n = 0; n<nodes.size(); n++)
  {
    if (nodes[n].found == 1)
      continue;
    bool dominated = positive ? root.somethingLarger(nodes, nodes[n].pos) : root.somethingSmaller(nodes, nodes[n].pos);
    if (dominated)
      nodes[n].found = 1;
    else
      front.push_back(nodes[n].pos);
  }
}

void Terrain::extract(const Cloud &cloud, const std::string &file_prefix, double width, bool verbose)
{
  // based on: Algorithms and Analyses for Maximal Vector Computation. Godfrey
  // but modified to use an Octal Space Partition tree. (like a BSP tree, but divided into 8 axis aligned per node)

  const double root_half = sqrt(0.5);
  const double root_3 = sqrt(3.0);
  const double root_2 = sqrt(2.0);
  Eigen::Matrix3d mat;
  mat.row(0) = Eigen::Vector3d(root_2/root_3,0,1.0/root_3);
  mat.row(1) = Eigen::Vector3d(-root_half/root_3, -root_half, 1.0/root_3);
  mat.row(2) = Eigen::Vector3d(-root_half/root_3, root_half, 1.0/root_3);

  Eigen::Matrix3d imat = mat.inverse();
  std::vector<Vector4d> points;
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
//    if (!cloud.rayBounded(i))
//      continue;
    Eigen::Vector3d pos = mat * cloud.ends[i];
    points.push_back(Vector4d(pos[0], pos[1], pos[2], (double)points.size() + 0.5));
  }

  RGBA white;
  white.red = white.green = white.blue = white.alpha = 255;  

  std::vector<Eigen::Vector3d> median_surface;
  std::vector<Vector4d> pos, neg;

  // first find the lower bound
  getParetoFront(points, neg, false);

  // if there is no width then we are already finished...
  if (width == 0.0)
  {
    // then convert it into a mesh
    std::vector<Eigen::Vector3d> vecs(neg.size());
    std::vector<Eigen::Vector3d> vecs_flat(neg.size());
    for (size_t i = 0; i<neg.size(); i++)
    {
      vecs[i] = imat * Eigen::Vector3d(neg[i][0], neg[i][1], neg[i][2]);
      vecs_flat[i] = vecs[i];
      vecs_flat[i][2] = 0.0;
    }
    ConvexHull hull(vecs_flat);
    hull.growUpwards(0.01); // same as a Delauney triangulation
    hull.mesh().vertices() = vecs;
    writePlyMesh(file_prefix + "_terrain_mesh.ply", hull.mesh(), true);
    if (verbose)
    {
      Cloud local_cloud;
      double t = 0.0;
      for (auto &p: vecs)
        local_cloud.addRay(p, p, t++, white);
      local_cloud.save(file_prefix + "_terrain.ply");
    }
    return;
  }

  cropPointsAbove(points, neg, width);

  int iteration = 0;
  while ((points.size() + pos.size() + neg.size()) > 0) // hopefully only do this O(n^1/3) times
  {
    std::vector<Vector4d> points_pos = points; 
    for (auto &p: neg)
      points_pos.push_back(Vector4d(p[0], p[1], p[2], (double)points_pos.size() + 0.5));
    pos.clear();
    if (iteration > 0)
    {
      std::vector<Vector4d> points_neg = points; 
      for (auto &p: pos)
        points_neg.push_back(Vector4d(p[0], p[1], p[2], (double)points_neg.size() + 0.5));
      neg.clear();
      getParetoFront(points_neg, neg, false);
    }
    getParetoFront(points_pos, pos, true); 
    std::vector<int> filled(points.size());
    memset(&filled[0], 0, sizeof(int)*filled.size());
    for (int i = (int)pos.size()-1; i>=0; i--) // auto &p: pos)
    {
      Vector4d p = pos[i];
      int index = (int)p[3];
      if (index >= (int)points.size()) // upwards has picked a downwards layer
      {
        median_surface.push_back(Eigen::Vector3d(p[0], p[1], p[2]));
        pos[i] = pos.back(); pos.pop_back();
      }
      else 
        filled[(int)p[3]] = 1;
    }
    for (int i = (int)neg.size()-1; i>=0; i--) // auto &n: neg)
    {
      Vector4d n = neg[i];
      int index = (int)n[3];
      if (index >= (int)points.size()) // downwards has picked an upwards layer
      {
   //     median_surface.push_back(Eigen::Vector3d(n[0], n[1], n[2]));
        neg[i] = neg.back(); neg.pop_back();
      }
      else
      {
        if (filled[index])
        {
          median_surface.push_back(Eigen::Vector3d(n[0], n[1], n[2]));
          neg[i] = neg.back(); neg.pop_back();        
        }
        filled[index] = 2;
      }
    }
    for (int i = (int)pos.size()-1; i>=0; i--) // auto &p: pos)
    {
      int index = (int)pos[i][3];
      if (filled[index] == 2)
      {
        pos[i] = pos.back(); pos.pop_back();
      }
    }
    std::cout << " points size: " << points.size() << " pos size: " << pos.size() << ", neg size: " << neg.size() << ", output size: " << median_surface.size() << std::endl;
    // now remove all the filled points
    for (int i = (int)points.size()-1; i>=0; i--) // could do this a little quicker
    {
      if (filled[i])
      {
        points[i] = points.back(); points.pop_back();
        points[i][3] = (double)i + 0.5; // readjust the indexing
      }
    }
    iteration++;
  }
  {
    std::vector<Eigen::Vector3d> vecs(median_surface.size());
    std::vector<Eigen::Vector3d> vecs_flat(median_surface.size());
    for (size_t i = 0; i<median_surface.size(); i++)
    {
      vecs[i] = imat * median_surface[i];
      vecs_flat[i] = vecs[i];
      vecs_flat[i][2] = 0.0;
    }
    ConvexHull hull(vecs_flat);
    hull.growUpwards(0.01); // same as a Delauney triangulation
    hull.mesh().vertices() = vecs;
    writePlyMesh(file_prefix + "_terrain_mesh.ply", hull.mesh(), true);
    if (verbose)
    {
      Cloud local_cloud;
      double t = 0.0;
      for (auto &p: vecs)
        local_cloud.addRay(p, p, t++, white);
      local_cloud.save(file_prefix + "_terrain.ply");
    }
  }
}

} // ray
