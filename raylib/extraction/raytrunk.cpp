// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrunk.h"
#include <nabo/nabo.h>
#include <map>
#include <queue>
#include "../raycuboid.h"
#include "../raygrid.h"

namespace ray
{

Trunk::Trunk()
  : centre(0, 0, 0)
  , radius(0)
  , score(0)
  , last_score(0)
  , length(0)
  , actual_length(0)
  , ground_height(0)
  , dir(0, 0, 0)
  , parent(-1)
  , active(true)
{}

// return the points that overlap this trunk, using the grid as an acceleration structure
std::vector<Eigen::Vector3d> Trunk::getOverlappingPoints(const Grid<Eigen::Vector3d> &grid, double spacing)
{
  std::vector<Eigen::Vector3d> points;
  // get grid bounds
  const Eigen::Vector3d base = centre - 0.5 * length * dir;
  const Eigen::Vector3d top = centre + 0.5 * length * dir;
  const double outer_radius = (radius + spacing) * boundary_radius_scale;
  const Eigen::Vector3d rad(outer_radius, outer_radius, outer_radius);
  const Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);
  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  mins = maxVector(mins, Eigen::Vector3i(0, 0, 0));
  const Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1, 1, 1);
  maxs = minVector(maxs, min_dims);

  // iterate over the cells in the bounds
  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0] <= maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1] <= maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2] <= maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        // iterate over the points in the cell
        for (auto &pos : cell.data)
        {
          // intersect against the trunk cylinder
          Eigen::Vector3d p = pos - centre;
          const double h = p.dot(dir);
          if (std::abs(h) > length * 0.5)
          {
            continue;
          }
          p -= dir * h;
          const double dist2 = p.squaredNorm();
          if (dist2 <= outer_radius * outer_radius)
          {
            points.push_back(pos);
          }
        }
      }
    }
  }
  return points;
}

// estimate the pose (centre and direction) of the trunk, from the set of points
void Trunk::estimatePose(const std::vector<Eigen::Vector3d> &points)
{
  centre = mean(points);
  // get teh scatter matrix from the points
  Eigen::Matrix3d scatter;
  scatter.setZero();
  for (const auto &point : points)
  {
    scatter += (point - centre) * (point - centre).transpose();
  }
  scatter /= static_cast<double>(points.size());
  // calculate an eigendecomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
  // the eigenvector with the largest eigenvalue (the long direction of the ellipsoid)
  // is the chosed trunk direction
  dir = eigen_solver.eigenvectors().col(2);
}

// improve the trunk's direction (dir) vector using the nearby set of points
void Trunk::updateDirection(const std::vector<Eigen::Vector3d> &points)
{
  // structure used for least squares best fit line
  struct Accumulator
  {
    Accumulator()
      : weight(0)
      , x(0)
      , abs_x(0)
      , x2(0)
      , y(0, 0)
      , xy(0, 0)
    {}
    double weight;
    double x;
    double abs_x;
    double x2;
    Eigen::Vector2d y;
    Eigen::Vector2d xy;
  };
  // start by taking two orthogonal axes to the current trunk dir
  const Eigen::Vector3d ax1 = Eigen::Vector3d(1, 2, 3).cross(dir).normalized();
  const Eigen::Vector3d ax2 = ax1.cross(dir);
  Accumulator sum;
  // for each point
  for (size_t i = 0; i < points.size(); i++)
  {
    // work relative to the current trunk centre
    const Eigen::Vector3d to_point = points[i] - centre;
    Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));

    // remove radius. If radius_removal_factor=0 then half-sided trees will have estimated trunk centred on that edge
    //                If radius_removal_factor=1 then v thin trunks may accidentally get a radius and it won't shrink
    //                down
    const double radius_removal_factor = 0.5;
    offset -= offset * radius_removal_factor * radius / offset.norm();

    const double h = to_point.dot(dir);
    sum.x += h;                // shift laong trunk length
    sum.y += offset;           // lateral offset relative to trunk
    sum.xy += h * offset;      // to obtain lean relative to trunk
    sum.x2 += h * h;           // to obtain spread of points along trunk
    sum.abs_x += std::abs(h);  // to obtain length of points
    sum.weight++;
  }
  const double n = sum.weight;
  centre += dir * (sum.x / n);

  // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
  Eigen::Vector2d sXY = sum.xy - sum.x * sum.y / n;
  const double sXX = sum.x2 - sum.x * sum.x / n;
  const double eps = 1e-10; // minimal value to use for the division
  if (std::abs(sXX) > eps)
  {
    sXY /= sXX;
  }

  dir = (dir + ax1 * sXY[0] + ax2 * sXY[1]).normalized();
  actual_length = 4.0 * (sum.abs_x / n);
  length = std::max(actual_length, 2.0 * radius * trunk_height_to_width);  // avoid getting too short for its width
}

// update estimation of the centre of the cylinder's circle
// taking a mean of the points is insufficient to estimate the centre when trunks are scanned
// from a single side. Instead we project the points to a paraboloid and look for a plane of
// best fit through this paraboloid.
void Trunk::updateCentre(const std::vector<Eigen::Vector3d> &points)
{
  // obtain two orthogonal vectors to the trunk's length
  const Eigen::Vector3d ax1 = Eigen::Vector3d(1, 2, 3).cross(dir).normalized();
  const Eigen::Vector3d ax2 = ax1.cross(dir);

  std::vector<Eigen::Vector3d> ps(points.size());
  Eigen::Vector3d mean_p(0, 0, 0);
  // for each point
  for (size_t i = 0; i < points.size(); i++)
  {
    // get offset relative to current trunk estimate
    const Eigen::Vector3d to_point = points[i] - centre;
    const Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));
    // make it relative to the trunk radius too
    const Eigen::Vector2d xy = offset / radius;
    const double l2 = xy.squaredNorm();
    // project the points to a paraboloid that has gradient 1 at 1
    const Eigen::Vector3d point(xy[0], xy[1], 0.5 * l2);
    ps[i] = point;
    mean_p += point;
  }
  mean_p /= static_cast<double>(points.size());
  // an accumulation structure for planes of best fit
  struct Acc
  {
    Acc() { x2 = y2 = xy = xz = yz = 0; }
    double x2, y2, xy, xz, yz;
  };
  Acc plane;
  for (auto &p : ps)
  {
    // accumulate the plane of best fit's parameters
    Eigen::Vector3d q = p - mean_p;
    plane.x2 += q[0] * q[0];
    plane.y2 += q[1] * q[1];
    plane.xy += q[0] * q[1];
    plane.xz += q[0] * q[2];
    plane.yz += q[1] * q[2];
  }
  // calculate the plane of best fit as local A,B values
  const double A = (plane.xz * plane.y2 - plane.yz * plane.xy) / (plane.x2 * plane.y2 - plane.xy * plane.xy);
  const double B = (plane.yz - A * plane.xy) / plane.y2;
  Eigen::Vector2d shift(A, B);
  const double shift2 = shift.squaredNorm();
  if (shift2 > 1.0)  // don't shift more than one radius each iteration
  {
    shift /= std::sqrt(shift2);
  }

  // convert this plane back into a lateral shift in the location of the trunk
  centre += (ax1 * shift[0] + ax2 * shift[1]) * radius;
}

// estimate the radius of the trunk from the points
void Trunk::updateRadius(const std::vector<Eigen::Vector3d> &points)
{
  // radius is just the mean radius relative to the previously calculated trunk pose
  double rad = 0;
  for (const auto &point : points)
  {
    const Eigen::Vector3d to_point = point - centre;
    rad += (to_point - dir * dir.dot(to_point)).norm();
  }
  const double n = static_cast<double>(points.size());
  radius = rad / n;
}


// calculate a score to represent how well the trunk fits to the points
void Trunk::updateScore(const std::vector<Eigen::Vector3d> &points)
{
  std::vector<double> scores(points.size());
  // now calculate the mean difference of the points from the expected trunk surface
  double raddiff = 0.0;
  for (const auto &point : points)
  {
    const Eigen::Vector3d to_point = point - centre;
    double rad = (to_point - dir * dir.dot(to_point)).norm();
    raddiff += std::abs(rad - radius);
  }
  // we subtract 4 here because this ia sample mean where 4 points were already used
  // to estimate the centre and lean of the trunk
  const double num_points = static_cast<double>(points.size()) - 4.0;
  const double variation = raddiff / num_points;  // end part gives sample variation

  last_score = score;
  // the variation per trunk length is scale invariant.
  // We use the reciprocal for the score
  score = actual_length / variation;
}

}  // namespace ray
