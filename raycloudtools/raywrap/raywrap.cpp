// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/rayconcavehull.h"
#include "raylib/rayconvexhull.h"
#include "raylib/raydebugdraw.h"
#include "raylib/rayply.h"
#include "raylib/rayutils.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>

// FIXME: Windows compatibility
#include <getopt.h>

void usage(int exitCode = 0)
{
  std::cout << "Extracts the ground surface as a mesh." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raywrap raycloud upwards 1.0 - wraps raycloud from the bottom upwards, or: downwards, inwards, outwards"
       << std::endl;
  std::cout << "                               the 1.0 is the maximum curvature to bend to" << std::endl;
  std::cout << "--full                       - the full (slower) method accounts for overhangs." << std::endl;
  exit(exitCode);
}

int main(int argc, char *argv[])
{
  if (argc < 4 || argc > 5)
    usage();

  ray::DebugDraw::init(argc, argv, "ConcaveHull");
  std::string file = argv[1];
  std::string type = argv[2];
  double maximum_curvature = std::stod(argv[3]);
  bool overhangs = false;
  if (argc == 5)
  {
    if (std::string(argv[4]) == "--full" or std::string(argv[4]) == "-f")
      overhangs = true;
    else
      usage();
  }

  ray::Cloud cloud;
  cloud.load(file);
  cloud.removeUnboundedRays();
  if (file.substr(file.length() - 4) == ".ply")
    file = file.substr(0, file.length() - 4);

  if (overhangs)
  {
    ray::ConcaveHull concave_hull(cloud.ends);
    if (type == "inwards")
      concave_hull.growInwards(maximum_curvature);
    else if (type == "outwards")
      concave_hull.growOutwards(maximum_curvature);
    else if (type == "upwards")
      concave_hull.growUpwards(maximum_curvature);
    else if (type == "downwards")
      concave_hull.growTopDown(maximum_curvature);
    else
      usage();

    ray::Mesh mesh;
    mesh.vertices = concave_hull.vertices;
    int num_bads = 0;
    for (auto &face : concave_hull.surface)
    {
      Eigen::Vector3d centroid(0, 0, 0);
      ConcaveHull::Tetrahedron &tetra = concave_hull.tetrahedra[face.tetrahedron];
      Eigen::Vector3i tri_verts = concave_hull.triangles[face.triangle].vertices;
      if (tri_verts[0] == -1)
        std::cout << "bad vertices in the surface" << std::endl;
      if (tetra.vertices[0] != -1)
      {
        for (int i = 0; i < 4; i++) centroid += concave_hull.vertices[tetra.vertices[i]] / 4.0;
        Eigen::Vector3d vs[3];
        for (int i = 0; i < 3; i++) vs[i] = concave_hull.vertices[tri_verts[i]];
        Eigen::Vector3d normal = (vs[2] - vs[0]).cross(vs[1] - vs[0]);
        if ((centroid - vs[0]).dot(normal) < 0.0)
          swap(tri_verts[1], tri_verts[2]);
      }
      else
        num_bads++;
      mesh.index_list.push_back(tri_verts);
    }
    if (num_bads > 0)
      std::cout << "number of surfaces that didn't have enough information to orient: " << num_bads << std::endl;
    writePlyMesh(file + "_mesh.ply", mesh, true);
  }
  else
  {
    ray::ConvexHull convexHull(cloud.ends);
    if (type == "inwards")
      convexHull.growInwards(maximum_curvature);
    else if (type == "outwards")
      convexHull.growOutwards(maximum_curvature);
    else if (type == "upwards")
      convexHull.growUpwards(maximum_curvature);
    else if (type == "downwards")
      convexHull.growTopDown(maximum_curvature);
    else
      usage();

    writePlyMesh(file + "_mesh.ply", convexHull.mesh, true);
  }


  std::cout << "Completed, output: " << file << "_mesh.ply" << std::endl;
  return 1;
}
