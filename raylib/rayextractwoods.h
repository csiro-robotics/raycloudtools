#pragma once
#include "standardIncludes.h"
#include "tools/raycloud.h"

#define LENGTH
struct Trunk
{
  Vector3d centre; // height is midway up trunk
  double radius;
  double score;
  double weight;
  double thickness;
  double length; 
  Vector2d lean;
};

struct Ray
{
  Ray(const Vector3d &start, const Vector3d &pos) : start(start), pos(pos) {}
  Vector3d start, pos;
};

struct Cell
{
  vector<Ray> rays;
  Vector2d minBound;
  double height;
};

struct Wood
{
  Wood(const RayCloud &cloud, double midRadius, double heightRange, const vector<Vector2d> &markedPoints);
  vector<Cell> grid;
  vector<Trunk> trunks;
  double width;
  Vector2i minBound, size;
};

struct Accumulator
{
  Accumulator(): weight(0), x(0), y(0,0), xy(0,0), x2(0), radius(0), radius2(0), z(0,0), xz(0,0) {}
  Accumulator &operator =(const Accumulator &a)
  {
    weight = a.weight;
    x = a.x;
    y = a.y;
    xy = a.xy;
    x2 = a.x2;
    radius = a.radius;
    radius2 = a.radius2;
    z = a.z;
    xz = a.xz;
    return *this;
  }
  Accumulator operator -(const Accumulator &a)
  {
    Accumulator res;
    res.weight = weight - a.weight;
    res.x = x - a.x;
    res.y = y - a.y;
    res.xy = xy - a.xy;
    res.x2 = x2 - a.x2;
    res.radius = radius - a.radius;
    res.radius2 = radius2 - a.radius2;
    res.z = z - a.z;
    res.xz = xz - a.xz;
    return res;
  }
  double weight;
  double x;
  Vector2d y;
  Vector2d xy;
  double x2;
  double radius;
  double radius2;
  Vector2d z;
  Vector2d xz;
};
