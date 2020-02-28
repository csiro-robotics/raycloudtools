#pragma once
#include "rayutils.h"
#include <set>

namespace RAY
{
struct ConvexHull
{
  ConvexHull(const vector<Vector3d> &points);

  void growOutwards(double maxCurvature);
  void growInwards(double maxCurvature);
  void growInDirection(double maxCurvature, const Vector3d &dir);
  void growUpwards(double maxCurvature){ growInDirection(maxCurvature, Vector3d(0,0,1)); }
  void growTopDown(double maxCurvature){ growInDirection(maxCurvature, Vector3d(0,0,-1)); }

  vector<Vector3d> points;
  vector<Vector3d> vertices;
protected:
  void ConvexHull::construct(const vector<Vector3d> &points);
};
}