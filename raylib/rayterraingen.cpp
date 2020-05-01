// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterraingen.h"
using namespace ray;
using namespace std;
using namespace Eigen;

struct PlanarWave
{
  PlanarWave() {}
  PlanarWave(const Vector2d &dir, double amplitude)
    : dir(dir)
    , amplitude(amplitude)
  {}
  Vector2d dir;  // includes wave frequency
  double amplitude;
};

// Some outdoor hilly terrain
void TerrainGen::generate()
{
  // 1. build the surface
  vector<PlanarWave> waves;
  for (int i = 0; i < 10; i++)
  {
    double amplitude = pow(0.8, i);
    Vector2d dir(random(-1.0, 1.0), random(-1.0, 1.0));
    dir.normalize();
    PlanarWave wave(0.1 * dir / amplitude, amplitude);
    waves.push_back(wave);
  }

  for (double time = 0.0; time < 100.0; time += 0.01)
  {
    // Now generate a trajectory around the landscape
    double phase = time * 2.0 * kPi / 100.0;
    Vector2d traj_centre(random(-10.0, 10.0), random(-10.0, 10.0));
    double traj_radius = random(5.0, 7.5);
    double rad_scale = random(0.05, 0.3);
    double traj_radius2 = traj_radius * rad_scale;
    Vector2d traj_pos = traj_centre + traj_radius * Vector2d(sin(phase), cos(phase)) +
                        traj_radius2 * Vector2d(sin(phase / rad_scale), cos(phase / rad_scale));

    double floor_y = 0.0;
    for (auto &wave : waves) floor_y += wave.amplitude * sin(traj_pos.dot(wave.dir));

    double scan_height = 1.3;
    Vector3d start(traj_pos[0], traj_pos[1], floor_y + scan_height);
    Vector3d dir(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, -0.6));
    dir.normalize();

    // now project the ray onto the terrain... how?
    // well we use an iterative scheme where we add the distance to the ground iteratively as a next guess
    double range = scan_height;
    for (int i = 0; i < 10; i++)
    {
      Vector3d pos = start + range * dir;
      Vector2d p(pos[0], pos[1]);
      double floor_y = 0.0;
      for (auto &wave : waves) floor_y += wave.amplitude * sin(p.dot(wave.dir));
      range += pos[2] - floor_y;
    }

    const double range_noise = 0.03;
    Vector3d end = start + (range + random(-range_noise, range_noise)) * dir;
    ray_starts.push_back(start);
    ray_ends.push_back(end);

    // what I really need to do is to expose a distance function... then the ray tracer can just pick the minimum (when
    // it is a union)
  }
}
