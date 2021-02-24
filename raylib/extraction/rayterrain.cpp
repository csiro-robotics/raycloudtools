// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterrain.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayconvexhull.h"
#include "../rayprogress.h"
#include "../rayprogressthread.h"

#if RAYLIB_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#endif  // RAYLIB_WITH_TBB


namespace ray
{
struct Node
{
  Node()
  { 
    pos.setZero();
    found = 0;
    is_set = 0;
    for (int i = 0; i<2; i++)
      for (int j = 0; j<2; j++)
        for (int k = 0; k<2; k++)
          dir_ids[i][j][k] = -1;
  }
  Vector4d pos;
  int found;
  int is_set;
  int dir_ids[2][2][2];
 
  bool somethingSmaller(std::vector<Node> &nodes, const Vector4d &corner)
  {
    #define CONE_CHECK
    Vector4d dif = corner - pos;
    if (dif == Vector4d(0,0,0,0))
    {
    #if defined CONE_CHECK
      if (dir_ids[0][0][0] == -1)
      {
        return false;
      }
      return nodes[dir_ids[0][0][0]].somethingSmaller(nodes, corner);
    #else
      return dir_ids[0][0][0] != -1;
    #endif
    }
    int i = (int)(dif[0] > 0.0);
    int j = (int)(dif[1] > 0.0);
    int k = (int)(dif[2] > 0.0);
    #if defined CONE_CHECK
    static const double root_third = std::sqrt(1.0/3.0); 
    static const double cos_ang = std::sqrt(2.0/3.0);
    static const Eigen::Vector3d diagonal(root_third, root_third, root_third);
    #endif
    if (i==0 && j==0 && k==0) // corner is smaller, so deactivate current node
    {
      #if defined CONE_CHECK
      Eigen::Vector3d dir = -Eigen::Vector3d(dif[0],dif[1], dif[2]).normalized();
      if (dir.dot(diagonal) > cos_ang)
        found = 1;
      #else
      found = 1; 
      #endif
    }
    else if (i==1 && j==1 && k==1) // corner is larger, so this node is indeed smaller
    {
      #if defined CONE_CHECK
      Eigen::Vector3d dir = Eigen::Vector3d(dif[0],dif[1], dif[2]).normalized();
      if (dir.dot(diagonal) > cos_ang)
      {
        return true;
      }
      #else
      return true; // needs a cone check
      #endif
    }
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
        break;
      }
      head = new_head;
    }
  }
}

void Terrain::getParetoFront(const std::vector<Vector4d> &points, std::vector<Vector4d> &front)
{
  std::vector<Node> nodes;
  constructOctalSpacePartition(nodes, points);
  Node &root = nodes[0];

  Progress progress;
  ProgressThread progress_thread(progress);
  progress.begin("rays processed:", nodes.size());
  
  const auto process_rays = [&nodes, &root, &front, &progress](size_t n)
  {
    progress.increment();
    if (nodes[n].found == 1)
    {
      return;
    }
    bool dominated = root.somethingSmaller(nodes, nodes[n].pos);
    if (dominated)
      nodes[n].found = 1;
    else
      nodes[n].is_set = 1;
  };
#if RAYLIB_WITH_TBB
  tbb::parallel_for<size_t>(0, nodes.size(), process_rays);
#else  
  for (size_t n = 0; n<nodes.size(); n++)
  {
    process_rays(n);
  }
#endif
  for (auto &node: nodes)
  {
    if (node.is_set)
    {
      front.push_back(node.pos);    
    }
  }
  progress.end();
  progress_thread.requestQuit();
  progress_thread.join();  
}

void Terrain::extract(const Cloud &cloud, const std::string &file_prefix, double gradient, bool verbose)
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
  double count = 0.5;
  double grad_scale = gradient / std::sqrt(2.0);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3d p = cloud.ends[i];
    p[2] /= grad_scale;
    Eigen::Vector3d pos = mat * p;
    points.push_back(Vector4d(pos[0], pos[1], pos[2], count));
    count++;
  }


  std::vector<Eigen::Vector3d> median_surface;
  std::vector<Vector4d> neg;

  // first find the lower bound
  getParetoFront(points, neg);
  std::cout << "number of wrap points: " << neg.size() << std::endl;

  // if there is no width then we are already finished...
  // then convert it into a mesh
  std::vector<Eigen::Vector3d> vecs(neg.size());
  std::vector<Eigen::Vector3d> vecs_flat(neg.size());
  for (size_t i = 0; i<neg.size(); i++)
  {
    vecs[i] = imat * Eigen::Vector3d(neg[i][0], neg[i][1], neg[i][2]);
    vecs[i][2] *= grad_scale;
    vecs_flat[i] = vecs[i];
    vecs_flat[i][2] = 0.0;
  }
  ConvexHull hull(vecs_flat);
  hull.growUpwards(0.01); // same as a Delauney triangulation
  hull.mesh().vertices() = vecs;
  writePlyMesh(file_prefix + "_mesh.ply", hull.mesh(), true);
  if (verbose)
  {
    RGBA white;
    white.red = white.green = white.blue = white.alpha = 255;  
    Cloud local_cloud;
    double t = 0.0;
    for (auto &p: vecs)
      local_cloud.addRay(p, p, t++, white);
    local_cloud.save(file_prefix + "_terrain.ply");
  }
}
} // ray
