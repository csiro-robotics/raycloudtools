// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayterraingen.h"
using namespace RAY;
using namespace std;
using namespace Eigen;

struct PlanarWave
{
  PlanarWave() {}
  PlanarWave(const Vector2d &dir, double amplitude) : dir(dir), amplitude(amplitude) {}
  Vector2d dir; // includes wave frequency
  double amplitude;
};

// Some outdoor hilly terrain
void TerrainGen::generate()
{
  // 1. build the surface
  vector<PlanarWave> waves;
  for (int i = 0; i<10; i++)
  {
    double amplitude = pow(0.8, i);
    Vector2d dir(random(-1.0, 1.0), random(-1.0, 1.0));
    dir.normalize();
    PlanarWave wave(0.1*dir/amplitude, amplitude);
    waves.push_back(wave);
  }

  for (double time = 0.0; time < 100.0; time+=0.01)
  {
    // Now generate a trajectory around the landscape
    double phase = time * 2.0*pi / 100.0;
    Vector2d trajCentre(random(-10.0, 10.0), random(-10.0, 10.0));
    double trajRadius = random(5.0, 7.5);
    double radScale = random(0.05, 0.3);
    double trajRadius2 = trajRadius * radScale;
    Vector2d trajPos = trajCentre + trajRadius*Vector2d(sin(phase), cos(phase)) + trajRadius2*Vector2d(sin(phase/radScale), cos(phase/radScale));
    
    double floorY = 0.0;
    for (auto &wave: waves)
      floorY += wave.amplitude * sin(trajPos.dot(wave.dir)); 
    
    double scanHeight = 1.3;
    Vector3d start(trajPos[0], trajPos[1], floorY + scanHeight);
    Vector3d dir(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, -0.6));
    dir.normalize();

    // now project the ray onto the terrain... how?
    // well we use an iterative scheme where we add the distance to the ground iteratively as a next guess
    double range = scanHeight;
    for (int i = 0; i<10; i++)
    {
      Vector3d pos = start + range*dir;
      Vector2d p(pos[0], pos[1]);
      double floorY = 0.0;
      for (auto &wave: waves)
        floorY += wave.amplitude * sin(p.dot(wave.dir)); 
      range += pos[2] - floorY;
    }

    const double rangeNoise = 0.03;
    Vector3d end = start + (range + random(-rangeNoise, rangeNoise))*dir;
    rayStarts.push_back(start);
    rayEnds.push_back(end);

    // what I really need to do is to expose a distance function... then the ray tracer can just pick the minimum (when it is a union)
  }

}
