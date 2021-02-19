// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterraingen.h"
#include "raymesh.h"
#include "rayply.h"

namespace ray
{
class PlanarWave
{
public:
  PlanarWave() {}
  PlanarWave(const Eigen::Vector2d &dir, double amplitude)
    : dir(dir)
    , amplitude(amplitude)
  {}
  Eigen::Vector2d dir;  // includes wave frequency
  double amplitude;
};

static const double range_noise = 0.03;

// Some outdoor hilly terrain
void TerrainGen::generate()
{
  // 1. build the surface
  std::vector<PlanarWave> waves;
  for (int i = 0; i < 10; i++)
  {
    double amplitude = pow(0.8, i);
    Eigen::Vector2d dir(random(-1.0, 1.0), random(-1.0, 1.0));
    dir.normalize();
    PlanarWave wave(0.1 * dir / amplitude, amplitude);
    waves.push_back(wave);
  }

  for (double time = 0.0; time < 100.0; time += 0.01)
  {
    // Now generate a trajectory around the landscape
    double phase = time * 2.0 * kPi / 100.0;
    Eigen::Vector2d traj_centre(random(-10.0, 10.0), random(-10.0, 10.0));
    double traj_radius = random(5.0, 7.5);
    double rad_scale = random(0.05, 0.3);
    double traj_radius2 = traj_radius * rad_scale;
    Eigen::Vector2d traj_pos = traj_centre + traj_radius * Eigen::Vector2d(sin(phase), cos(phase)) +
                        traj_radius2 * Eigen::Vector2d(sin(phase / rad_scale), cos(phase / rad_scale));

    double floor_y = 0.0;
    for (auto &wave : waves) floor_y += wave.amplitude * sin(traj_pos.dot(wave.dir));

    double scan_height = 1.3;
    Eigen::Vector3d start(traj_pos[0], traj_pos[1], floor_y + scan_height);
    Eigen::Vector3d dir(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, -0.6));
    dir.normalize();

    // now project the ray onto the terrain... how?
    // well we use an iterative scheme where we add the distance to the ground iteratively as a next guess
    double range = scan_height;
    for (int i = 0; i < 10; i++)
    {
      Eigen::Vector3d pos = start + range * dir;
      Eigen::Vector2d p(pos[0], pos[1]);
      double floor_y = 0.0;
      for (auto &wave : waves) floor_y += wave.amplitude * sin(p.dot(wave.dir));
      range += pos[2] - floor_y;
    }

    Eigen::Vector3d end = start + (range + random(-range_noise, range_noise)) * dir;
    ray_starts_.push_back(start);
    ray_ends_.push_back(end);

    // what I really need to do is to expose a distance function... then the ray tracer can just pick the minimum (when
    // it is a union)
  }
}

bool TerrainGen::generateFromFile(const std::string &filename)
{
  // simply rain down rays from a height, onto the mesh at a given density
  const double ray_height = 1.0; 
  const double ray_density = 400; // ray ends per square metre
  Mesh mesh;
  if (!readPlyMesh(filename, mesh))
  {
    return false;
  }
  std::vector<Eigen::Vector3i> &triangles = mesh.indexList();
  std::vector<Eigen::Vector3d> &vertices = mesh.vertices();
  // for each triangle, add rays in proportion to its horizontal area
  double area_covered = 0;
  double total_area = 0.0;
  const double area_per_ray = 1.0 / ray_density; 
  for (auto &triangle: triangles)
  {
    Eigen::Vector3d vs[3] = {vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]};
    vs[0][2] = vs[1][2] = vs[2][2] = 0.0; // flatten
    double area = std::abs(((vs[1] - vs[0]).cross(vs[2] - vs[0]))[2]) / 2.0;
    total_area += area;
    while (area_covered < total_area)
    {
      double ws[2] = {random(0,1), random(0,1)}; // weights give a random location
      if (ws[0] + ws[1] > 1.0)
      { 
        ws[0] = 1.0 - ws[0];
        ws[1] = 1.0 - ws[1];
      } 
      Eigen::Vector3d v0 = vertices[triangle[0]];
      Eigen::Vector3d v1 = vertices[triangle[1]];
      Eigen::Vector3d v2 = vertices[triangle[2]];
      Eigen::Vector3d ray_end = v0 + (v1-v0)*ws[0] + (v2-v0)*ws[1];

      ray_end[2] += random(-range_noise, range_noise);
      Eigen::Vector3d ray_syart = ray_end + Eigen::Vector3d(0,0,ray_height);
      ray_starts_.push_back(ray_syart);
      ray_ends_.push_back(ray_end);
      area_covered += area_per_ray;
    }
  }
  return true;
}

} // ray