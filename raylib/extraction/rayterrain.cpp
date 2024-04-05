// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterrain.h"
#include "../rayconvexhull.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayprogress.h"
#include "../rayprogressthread.h"

#if RAYLIB_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB
static int num_visits = 0;
static int num_cone_tests = 0;


namespace ray
{
/// The node structure used in calculating the pareto front
struct Node
{
  Node()
  {
    pos.setZero();
    found = 0;
    is_set = 0;
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        for (int k = 0; k < 2; k++)
        {
          dir_ids[i][j][k] = -1;
        }
      }
    }
  }
  Vector4d pos;
  int found;
  int is_set;
  int dir_ids[2][2][2];

  // returns whether there is a smaller node than the supplied @c corner point
  bool somethingSmaller(std::vector<Node> &nodes, const Vector4d &corner)
  {
    num_visits++;
// This checks in a cone rather than just the corner of a cube shape that you would get
// from a raw Pareto front calculation    
#define CONE_CHECK  
    Vector4d dif = corner - pos;
    if (dif == Vector4d(0, 0, 0, 0))
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
    const int i = static_cast<int>(dif[0] > 0.0);
    const int j = static_cast<int>(dif[1] > 0.0);
    const int k = static_cast<int>(dif[2] > 0.0);
#if defined CONE_CHECK
    static const double root_third = std::sqrt(1.0 / 3.0);
    static const double cos_ang = std::sqrt(2.0 / 3.0);
    static const Eigen::Vector3d diagonal(root_third, root_third, root_third);
#endif
    if (i == 0 && j == 0 && k == 0)  // corner is smaller, so deactivate current node
    {
#if defined CONE_CHECK
      num_cone_tests++;
      const Eigen::Vector3d dir = -Eigen::Vector3d(dif[0], dif[1], dif[2]).normalized();
      if (dir.dot(diagonal) > cos_ang)
        found = 1;
#else
      found = 1;
#endif
    }
    else if (i == 1 && j == 1 && k == 1)  // corner is larger, so this node is indeed smaller
    {
#if defined CONE_CHECK
      num_cone_tests++;
      const Eigen::Vector3d dir = Eigen::Vector3d(dif[0], dif[1], dif[2]).normalized();
      if (dir.dot(diagonal) > cos_ang)
      {
        return true;
      }
#else
      return true;  // needs a cone check
#endif
    }
    for (int I = 0; I <= i; I++)
    {
      for (int J = 0; J <= j; J++)
      {
        for (int K = 0; K <= k; K++)
        {
          if (dir_ids[I][J][K] != -1)
          {
            if (nodes[dir_ids[I][J][K]].somethingSmaller(nodes, corner))
            {
              return true;
            }
          }
        }
      }
    }
    return false;
  }
};

/// construct an octal space partition tree
void constructOctalSpacePartition(std::vector<Node> &nodes, std::vector<Vector4d> points)
{
  nodes.resize(points.size());
  int i = 0;
  while (points.size() > 0)
  {
    const int ind = rand() % static_cast<int>(points.size());
    nodes[i++].pos = points[ind];
    points[ind] = points.back();
    points.pop_back();
  }
  // n log n on average, for each node, trace the tree to add correct direction ids
  for (size_t n = 1; n < nodes.size(); n++)
  {
    const Eigen::Vector4d &pos = nodes[n].pos;
    int head = 0;
    for (;;)
    {
      const Vector4d dif = pos - nodes[head].pos;
      const int i = static_cast<int>(dif[0] > 0.0);
      const int j = static_cast<int>(dif[1] > 0.0);
      const int k = static_cast<int>(dif[2] > 0.0);
      const int new_head = nodes[head].dir_ids[i][j][k];
      if (new_head == -1)
      {
        nodes[head].dir_ids[i][j][k] = static_cast<int>(n);
        break;
      }
      head = new_head;
    }
  }
}

// get 3D pareto front, the 4D vectors' last element is its index, to aid with book keeping
void Terrain::getParetoFront(const std::vector<Vector4d> &points, std::vector<Vector4d> &front)
{
  // this is an acceleration structure for faster lookup
  std::vector<Node> nodes;
  constructOctalSpacePartition(nodes, points);
  Node &root = nodes[0];

  Progress progress;
  ProgressThread progress_thread(progress);
  progress.begin("rays processed:", nodes.size());

  const auto process_rays = [&nodes, &root, &front, &progress](size_t n) {
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
  #pragma omp parallel for
  for (size_t n = 0; n < nodes.size(); n++)
  {
    process_rays(n);
  }
#endif
  for (auto &node : nodes)
  {
    if (node.is_set)
    {
      front.push_back(node.pos);
    }
  }
  progress.end();
  progress_thread.requestQuit();
  progress_thread.join();
  std::cout << "number of rays: " << points.size() << ", number of visits: " << num_visits
            << ", number of cone tests: " << num_cone_tests << std::endl;
}

void Terrain::growUpwards(const std::vector<Eigen::Vector3d> &positions, double gradient)
{
#if RAYLIB_WITH_QHULL
  // The idea behind ground extraction is to tilt the upwards vector to the (1,1,1) direction then 
  // find the Pareto front in the three principle axes. https://en.wikipedia.org/wiki/Pareto_front
  //
  // Efficient Pareto front calculation is based on: Algorithms and Analyses for Maximal Vector Computation. Godfrey
  // but modified to use an Octal Space Partition tree. (like a BSP tree, but divided into 8 axis aligned per node)
  //
  // Parento fronts approximate a gradient constraint on the ground surface, but we modify the algorithm
  // in the CONE_CHECK code, to apply an exact gradient constraint.
  const double root_half = sqrt(0.5);
  const double root_3 = sqrt(3.0);
  const double root_2 = sqrt(2.0);

  // the following matrix is because we are rotating a vertical direction into the long-diagonal (1,1,1)
  // direction, because we are primarily treating the extraction problem as one of finding a
  // pareto front in 3D.
  Eigen::Matrix3d mat;
  mat.row(0) = Eigen::Vector3d(root_2 / root_3, 0, 1.0 / root_3);
  mat.row(1) = Eigen::Vector3d(-root_half / root_3, -root_half, 1.0 / root_3);
  mat.row(2) = Eigen::Vector3d(-root_half / root_3, root_half, 1.0 / root_3);
  Eigen::Matrix3d imat = mat.inverse();
  const double grad_scale = gradient / std::sqrt(2.0);

  // generate these transformed points as input
  std::vector<Eigen::Vector4d> points(positions.size());
  double count = 0.5;
  for (size_t i = 0; i < positions.size(); i++)
  {
    Eigen::Vector3d p = positions[i];
    p[2] /= grad_scale;
    const Eigen::Vector3d pos = mat * p;
    points[i] = Eigen::Vector4d(pos[0], pos[1], pos[2], count++);
  }

  std::vector<Vector4d> front;
  // then find pareto the lower bound of the points
  getParetoFront(points, front);
  std::cout << "number of pareto front points: " << front.size() << std::endl;

  // then convert it into a mesh
  std::vector<Eigen::Vector3d> vecs(front.size());
  std::vector<Eigen::Vector3d> vecs_flat(front.size());
  for (size_t i = 0; i < front.size(); i++)
  {
    // we convert the points back to world space
    vecs[i] = imat * Eigen::Vector3d(front[i][0], front[i][1], front[i][2]);
    vecs[i][2] *= grad_scale;
    // flattened points are used to get the Delauney triangulation below
    vecs_flat[i] = vecs[i];
    vecs_flat[i][2] = 0.0;
  }
  ConvexHull hull(vecs_flat);
  hull.growUpwards(0.01);  // same as a Delauney triangulation

  mesh_.indexList() = hull.mesh().indexList();
  mesh_.vertices() = vecs;
#endif
}

void Terrain::growDownwards(const std::vector<Eigen::Vector3d> &positions, double gradient)
{
  // as you might imagine, this is just the reverse of growupwards
  std::vector<Eigen::Vector3d> upsidedown_points = positions;
  for (auto &point : upsidedown_points)
  {
    point[2] = -point[2];
  }
  growUpwards(upsidedown_points, gradient);
  for (auto &point : mesh_.vertices())
  {
    point[2] = -point[2];
  }
}

// a faster version of the growupwards algorithm
void Terrain::growUpwardsFast(const std::vector<Eigen::Vector3d> &ends, double pixel_width,
                              const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double gradient)
{
#if RAYLIB_WITH_QHULL
  // the speed up is one of removing lots of 'above ground' points before running the growUpwards function
  // thereby making the problem size smaller.

  // the way we do this is to put the points into a 2D height grid, and do a course removal based on height

  // here we generate the grid and find the lowest points per cell
  const Eigen::Vector3d extent = max_bound - min_bound;
  const Eigen::Vector2i dims(static_cast<int>(std::ceil(extent[0] / pixel_width)),
                             static_cast<int>(std::ceil(extent[1] / pixel_width)));
  std::vector<Eigen::Vector3d> lowests(dims[0] * dims[1]);
  for (int i = 0; i < dims[0]; i++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      lowests[i + dims[0] * j] = Eigen::Vector3d(0, 0, std::numeric_limits<double>::max());
    }
  }
  for (size_t i = 0; i < ends.size(); i++)
  {
    const Eigen::Vector3d pos = (ends[i] - min_bound) / static_cast<double>(pixel_width);
    const int index = static_cast<int>(pos[0]) + dims[0] * static_cast<int>(pos[1]);
    if (ends[i][2] < lowests[index][2])
    {
      lowests[index] = ends[i];
    }
  }

  std::vector<Eigen::Vector3d> points;
  // then for each point
  for (size_t i = 0; i < ends.size(); i++)
  {
    const Eigen::Vector3d p = ends[i];
    // we get its cell index
    const Eigen::Vector3d point = (p - min_bound) / static_cast<double>(pixel_width);
    const int I = static_cast<int>(point[0]);
    const int J = static_cast<int>(point[1]);
    const int Imin = std::max(0, I - 1);
    const int Jmin = std::max(0, J - 1);
    const int Imax = std::min(I + 1, dims[0] - 1);
    const int Jmax = std::min(J + 1, dims[1] - 1);
    bool remove = false;
    // then within the Moore neighbourhood
    for (int x = Imin; x <= Imax; x++)
    {
      for (int y = Jmin; y <= Jmax; y++)
      {
        const int index = x + dims[0] * y;
        const Eigen::Vector3d d = p - lowests[index];
        const double dist2 = (d[0] * d[0] + d[1] * d[1]);
        // remove if neighbour points are lower
        if (d[2] > 0.0 && d[2] * d[2] > dist2)
        {
          remove = true;
          break;
        }
      }
      if (remove)
      {
        break;
      }
    }
    if (remove)
    {
      continue;
    }

    points.push_back(p);
  }
  std::cout << "size before: " << ends.size() << ", size after: " << points.size() << std::endl;

  growUpwards(points, gradient);
#endif
}

// Convert the @c cloud input to the mesh_ member variable. 
void Terrain::extract(const Cloud &cloud, const Eigen::Vector3d &offset, const std::string &file_prefix, double gradient, bool verbose)
{
#if RAYLIB_WITH_QHULL
  // preprocessing to make the cloud smaller.
  Eigen::Vector3d min_bound, max_bound;
  cloud.calcBounds(&min_bound, &max_bound);
  const double spacing = cloud.estimatePointSpacing();
  const double pixel_width = 2.0 * spacing;
  std::vector<Eigen::Vector3d> ends;
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (cloud.rayBounded(i))
    {
      ends.push_back(cloud.ends[i]);
    }
  }
  growUpwardsFast(ends, pixel_width, min_bound, max_bound, gradient);
  mesh_.reduce();  // remove disconnected vertices in the mesh
  mesh_.colours() = std::vector<RGBA>(mesh_.vertices().size(), RGBA::terrain());

  mesh_.translate(offset);
  writePlyMesh(file_prefix + "_mesh.ply", mesh_, true);
  if (verbose)  // debugging output
  {
    Cloud local_cloud;
    double t = 0.0;
    for (auto &p : mesh_.vertices())
    {
      local_cloud.addRay(p, p, t++, RGBA::white());
    }
    local_cloud.save(file_prefix + "_terrain.ply");
  }
#else
  std::cerr << "Error: extracting terrain requires QHull, see README instructions for installation" << std::endl;
#endif
}
}  // namespace ray
