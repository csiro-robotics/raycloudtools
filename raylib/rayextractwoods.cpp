// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayextractwoods.h"
#include "raydebugdraw.h"

namespace ray
{
#define LEAN

Eigen::Vector2d vector2d(const Eigen::Vector3d &v)
{
  return Eigen::Vector2d(v[0], v[1]);
}
Eigen::Vector3d vector3d(const Eigen::Vector2d &v, double z = 0)
{
  return Eigen::Vector3d(v[0], v[1], z);
}

#if OLD_WOODS
Wood::Wood(const Cloud &cloud, double midRadius, double heightRange, bool verbose)
{
  if (verbose)
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);

  Eigen::Vector3d min_bound = cloud.calcMinBound();
  Eigen::Vector3d max_bound = cloud.calcMaxBound();
  std::cout << "cloud from: " << min_bound.transpose() << " to: " << max_bound.transpose() << std::endl;
  
  const int minPointsPerTrunk = 15;
  const double minScore = 1.1;
  double maxRadius = midRadius * 2.0;

  width = (1.5 + 1e-5) * maxRadius; 
  Eigen::Vector2i minCell = Eigen::Vector2i(floor(min_bound[0] / width), floor(min_bound[1] / width));
  Eigen::Vector2i maxCell(floor(max_bound[0] / width), floor(max_bound[1] / width));
  Eigen::Vector2d minBound = Eigen::Vector2d(minCell[0], minCell[1]) * width; 
  size = Eigen::Vector2i(maxCell[0] + 1 - minCell[0], maxCell[1] + 1 - minCell[1]);
  grid.resize(size[0] * size[1]);
  for (int x = 0; x < size[0]; x++)
    for (int y = 0; y < size[1]; y++)
      grid[x + size[0] * y].minBound = minBound + Eigen::Vector2d(x, y) * width;

  // populate the cells
  int num_added = 0;
  for (int i = 0; i < (int)cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector2i index(floor((cloud.ends[i][0] - minBound[0]) / width), floor((cloud.ends[i][1] - minBound[1]) / width));
    Cell &cell = grid[index[0] + size[0] * index[1]];
    if (cell.rays.empty())
      num_added++;
    cell.rays.push_back(Ray(cloud.starts[i], cloud.ends[i]));
    cell.height += cloud.ends[i][2];
  }
  std::cout << "num cells added: " << num_added << std::endl;

  Trunk trunk;
  trunk.radius = midRadius;
  trunk.score = 0;
  for (auto &cell : grid)
  {
    cell.height /= (1e-10 + (double)cell.rays.size());
    if (cell.rays.size() >= minPointsPerTrunk)
    {
      trunk.centre = vector3d(cell.minBound + width * Eigen::Vector2d(0.5,0.5));
      trunk.length = heightRange / 3.0; // the catchment starts at twice this, so catches 2/3 of the passed in cloud
      trunk.centre[2] = cell.height;
      trunk.lean = Eigen::Vector2d(0,0);
      trunk.thickness = 0.0;
      trunks.push_back(trunk);
    }
  }
  std::cout << "number of trunks added initially: " << trunks.size() << std::endl;
  
 // DebugDraw::instance()->drawText(components(trunks, centre), vector<std::string>());

  
  const int maxIterations = 10;
  for (int it = 0; it < maxIterations; it++)
  {
    if (verbose)
      DebugDraw::instance()->drawTrunks(trunks);
    std::vector<Eigen::Vector3d> spheres;
    std::vector<Eigen::Vector3d> ringCentres;
    std::vector<Eigen::Vector3d> ringNormals;
    std::vector<double> ringRadii;
    double bestScore = 0.0;
    double radiusScale = 2.0 - 0.6*(double)it/(double)(maxIterations-1); // pull in the average, to get a better estimation after it converges
    double lengthScale = 2.0 - (double)it/(double)maxIterations;
    for (int tr = (int)trunks.size()-1; tr>=0; tr--)
    {
      Trunk &trunk = trunks[tr];
      // for each trunk candidate, find the (max) 4 cells it overlaps 
      double outerRadius = trunk.radius * radiusScale;
      Eigen::Vector2d minPos = vector2d(trunk.centre) - outerRadius*Eigen::Vector2d(1,1) - minBound;
      Eigen::Vector2d maxPos = vector2d(trunk.centre) + outerRadius*Eigen::Vector2d(1,1) - minBound;
      Eigen::Vector2i indexMin(std::max(0.0, floor(minPos[0] / width)), std::max(0.0, floor(minPos[1] / width)));
      Eigen::Vector2i indexMax(std::min((int)floor(maxPos[0] / width), size[0]-1), std::min((int)floor(maxPos[1] / width), size[1]-1));
      Eigen::Vector3d average(0, 0, 0);
      trunk.weight = 0;
      std::vector<Eigen::Vector3d> ps;
      const bool accountForInwardRays = it>0;
      std::vector<Ray> rays[5];

      double step = trunk.length/3.0;
      double heights[6] = {-lengthScale, -1.0, -1.0/lengthScale, 1.0/lengthScale, 1.0, lengthScale}; 
      for (int i = 0; i<6; i++)
        heights[i] = trunk.centre[2] + heights[i]*(0.5*trunk.length);
      double minh = 1e10;
      for (int x = indexMin[0]; x <= indexMax[0]; x++)
      {
        for (int y = indexMin[1]; y <= indexMax[1]; y++)
        {
          Cell &cell = grid[x + size[0] * y];
          minh = std::min(minh, cell.height);
          for (auto &rawRay : cell.rays)
          {
            Ray ray = rawRay;
#if defined(LEAN) // unlean time
            Eigen::Vector3d offset = vector3d(trunk.lean) * (ray.pos[2] - trunk.centre[2]); 
            ray.pos -= offset;
            ray.start -= offset; 
#endif
            Eigen::Vector2d xy = vector2d(ray.pos - trunk.centre)/trunk.radius;
            double l2 = xy.squaredNorm();
            if (l2 > sqr(outerRadius / trunk.radius))
              continue;
            if (ray.pos[2] > heights[0] && ray.pos[2] < heights[5])
            {
              // we're transforming around a Riemann sphere at trunk.centre with 90 degrees at trunk.radius
              Eigen::Vector3d point(2.0*xy[0] / (1.0 + l2), 2.0*xy[1] / (1.0 + l2), (1.0 - l2) / (1.0 + l2));
              // TODO: test using SVD instead of PCA (eigen decomp)
              ps.push_back(point);
              spheres.push_back((point + Eigen::Vector3d(0,0,0.5)) * trunk.radius + trunk.centre);
              
              average += point;
              trunk.weight++;
              int h = 0; 
              while (ray.pos[2] > heights[h+1])
                h++;
              if (h >= 5)
                std::cout << "huh?" << std::endl;
              rays[h].push_back(ray);
            }
          }
        }
      }
      
      if (trunk.weight == 0.0)
      {
        trunks[tr] = trunks.back();
        trunks.pop_back();
        continue;
      }
      average /= (double)trunk.weight;
      
      Eigen::Matrix3d scatter;
      scatter.setZero();
      for (auto &p : ps)
        scatter += (p-average) * (p-average).transpose();

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(scatter); 
      ASSERT(eigenSolver.info() == Success); // use EigenSolver if it isn't always self-adjoint
      Eigen::Vector3d eigenValues = eigenSolver.eigenvalues();
      Eigen::Matrix3d eigenVectors = eigenSolver.eigenvectors();

      Eigen::Vector3d normal = eigenVectors.col(0); // least eigenvaluestruct Ray
      if (normal.dot(average) < 0.0)
        normal = -normal;
      ringNormals.push_back(normal);
      Eigen::Vector2d flatNormal = vector2d(normal);

      // now work out centre and radius from projected stats
      double phi = std::atan2(flatNormal.norm(), normal[2]);
      double h = average.dot(normal);
      ringCentres.push_back((normal*h + Eigen::Vector3d(0,0,0.5)) * trunk.radius + trunk.centre);
      ringRadii.push_back(sqrt(1.0 - sqr(h)) * trunk.radius);
      double angle = std::acos(h);
      
      // make sure circle doesn't cross infinity (-z) line.
      double a = phi - angle;
      if (a < -kPi)
        a += 2.0*kPi;
      double b = phi + angle;
      if (b > kPi)
        b -= 2.0*kPi;
      
      double lMin = std::tan(std::min(a, b)/2.0) * trunk.radius;
      double lMax = std::tan(std::max(a, b)/2.0) * trunk.radius;
      double length = (lMax + lMin) / 2.0;
      flatNormal.normalize();
      double newRadius = (lMax - lMin) / 2.0;
      
      double intersectionArea;
      double radDiff = std::abs(newRadius - trunk.radius);
      // from https://www.xarg.org/2016/07/calculate-the-intersection-area-of-two-circles/
      if (abs(newRadius - trunk.radius) < 1e-5 && length < 1e-5)
        intersectionArea = kPi * sqr((newRadius + trunk.radius)/2.0);
      else if (length < radDiff)
        intersectionArea = kPi * sqr(std::min(newRadius, trunk.radius));
      else if (length > newRadius + trunk.radius)
        intersectionArea = 0;
      else
      {
        double R = trunk.radius;
        double r = newRadius;
        double d = length;
        double R2 = R*R;
        double r2 = r*r;
        double d2 = d*d;
        
        intersectionArea = r2*acos((d2+r2-R2)/(2*d*r))+R2*acos((d2+R2-r2)/(2*d*R)) - 
          0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
        
        if (intersectionArea > kPi * std::min(R2, r2) + 1e-9)
          std::cout << "weird intersection area" << std::endl;
      }
      double scoreScale = (kPi*sqr(trunk.radius) * kPi*sqr(newRadius))/(1e-10 + sqr(intersectionArea)); // this grows with the square of the distance between the disks
      const double pruningLateness = 1.0; // if you want it to prune more early and aggressively, try reducing this to 0.5 
      scoreScale = pow(scoreScale, pruningLateness);
      if (it == maxIterations-1)
        scoreScale = 1.0; // no advantage once the iterations are over
      
      trunk.centre += vector3d(flatNormal*length);
      trunk.radius = newRadius;
      
      std::vector<Ray> rays2;
     
      Accumulator cumulative[6];
      
      // now adjust the length and position
      for (int i = 1; i<6; i++)
      {
        cumulative[i] = cumulative[i-1];
        for (auto &ray: rays[i-1])
        {
          Eigen::Vector3d dif = ray.pos - trunk.centre;
          Eigen::Vector2d vec = vector2d(dif);
          double d = vec.norm();
          double dR = d - trunk.radius;
          Eigen::Vector2d pos = vec * dR/d;

          cumulative[i].y += pos;
          cumulative[i].x += dif[2];
          cumulative[i].xy += dif[2]*pos;
          cumulative[i].x2 += sqr(dif[2]);
          cumulative[i].z += -vec/d;
          cumulative[i].xz += -dif[2]*vec/d;
          cumulative[i].radius += dR;
          cumulative[i].radius2 += sqr(dR);
          cumulative[i].weight++;
        }
      }
      double maxScore = 0;
      int I = 1, J = 4;
      double heightVariance = 0;
      Eigen::Vector2d bestAbsLean(0,0);
      Eigen::Vector2d bestOffset(0,0);
      double bestRadius = 0;
      double bestThickness = 0;
      for (int i = 0; i<=3; i++)
      {
        for (int j = i+2; j<6; j++)
        {
          if ((double)(j-i)*step < 2.0*midRadius) // you can't make it smaller
            continue;
          // average variance inverted
          Accumulator sum = cumulative[j] - cumulative[i];
          if (sum.weight < minPointsPerTrunk)
            continue;
          double N = sum.weight;

          // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
          double sYY = sum.radius2 - sum.y.squaredNorm()/N;
          Eigen::Vector2d sXY = sum.xy - sum.x*sum.y/N;
          double sXX = sum.x2 - sum.x*sum.x/N;

          Eigen::Vector2d lean = sXY / sXX;
          
          Eigen::Vector2d absLean = lean + trunk.lean;
          Eigen::Vector2d oldLean = lean;
          double l = absLean.norm();
          double maxLean = 0.25;
          if (l > maxLean)
            absLean *= maxLean/l;          
          lean = absLean - trunk.lean;
          
          Eigen::Vector2d offset = sum.y/N - lean*sum.x/N;
          
          double thickness = sYY - lean.squaredNorm() * sXX;
          
#define PARALLEL_RADIUS // This adjusts the radius each segment, to get a smaller thickness. Otherwise the radius is adjusted after the segment is found
#if defined PARALLEL_RADIUS
          double dRadius = sum.radius/N + sum.z.dot(offset)/N + sum.xz.dot(lean)/N;
          thickness -= sqr(dRadius)*N;
          double radius = trunk.radius + dRadius; // we also quickly update the radius here, since length and lean may have changed
#endif      
          if (thickness < 0.0)
            std::cout << "weird, negative variance- " << oldLean.transpose() << std::endl;

          double score = N / (1e-10 + thickness);
          double hVariance = (sum.x2 / N) - sqr(sum.x / N);    
          score *= hVariance / (heights[j]-heights[i]); // this division ensures the scores are the same if the variance is proportional to the height of the trunk
          if (score > maxScore)
          {
            maxScore = score;
            heightVariance = hVariance;
            trunk.weight = N;
            bestAbsLean = absLean;
            bestOffset = offset;
#if defined PARALLEL_RADIUS
            bestRadius = radius;
            bestThickness = thickness;
#endif
            I=i;
            J=j;
          }
        }
      }
      trunk.centre += vector3d(bestOffset);
      trunk.lean = bestAbsLean;
      trunk.centre[2] = (heights[I] + heights[J])/2.0;
      trunk.length = heights[J]-heights[I];
      for (int i = I; i<J; i++)
        rays2.insert(rays2.end(), rays[i].begin(), rays[i].end());
      double thickness;
#if defined PARALLEL_RADIUS
      trunk.radius = bestRadius;
      thickness = bestThickness;
#endif
#define PERFECT_THICKNESS_AND_RADIUS // runs a last iteration on moved, tilted, lengthened cylinder, to get a more accurate thickness and radius
#if !defined PARALLEL_RADIUS || defined PERFECT_THICKNESS_AND_RADIUS
      double rad = 0;
      double rad2 = 0;
      for (auto &ray: rays2)
      {
        Eigen::Vector3d dif = ray.pos - trunk.centre;
        Eigen::Vector2d vec = vector2d(dif) - trunk.lean * dif[2];
        double d = vec.norm();
        rad += d - trunk.radius;
        rad2 += sqr(d - trunk.radius);
      }
      trunk.radius += rad/trunk.weight; // we also quickly update the radius here, since length and lean may have changed
      thickness = rad2 - sqr(rad) / trunk.weight; // this thickness is now relative to the new radius
#endif
      
      if (trunk.radius > midRadius * 3.0 || trunk.radius < midRadius / 3.0)
      {
        trunks[tr] = trunks.back();
        trunks.pop_back();
        continue;
      }
            
      if (accountForInwardRays)
      {
        for (auto &ray : rays2)
        {
          Eigen::Vector2d toStart = vector2d(ray.start - ray.pos);
          Eigen::Vector2d toCentre = vector2d(trunk.centre - ray.pos);
          if (toStart.dot(toCentre) > 0.0)
          {
            double dist1 = toCentre.norm();
            double dot = clamped(toCentre.dot(toStart) / toStart.squaredNorm(), 0.0, 1.0);
            double dist2 = (-toCentre + toStart * dot).norm();
            double delta = sqr(trunk.radius - dist2) - sqr(trunk.radius - dist1); // adjust the sum of square errors
            thickness += std::max(0.0, 4.0 * delta);
          }
        }
      }
      trunk.thickness = std::sqrt(thickness / trunk.weight);
      // TODO: try replacing max(eigen...) with sqr(radius) * trunk.weight. Should preference partial trees
      trunk.score = std::sqrt((heightVariance / trunk.length) * std::max(eigenValues[1], eigenValues[2]) / thickness); // max would do better with trees seen from one side, but more false positives
      // remove poorly performing trunks
      if (trunk.weight < minPointsPerTrunk || trunk.score*scoreScale < minScore)
      {
        trunks[tr] = trunks.back();
        trunks.pop_back();
        continue;
      }
      bestScore = std::max(bestScore, trunk.score);
    }
    std::cout << "best score: " << bestScore << std::endl;
    if (verbose)
    {
      DebugDraw::instance()->drawCloud(spheres, 1.0, 1);
      DebugDraw::instance()->drawRings(ringCentres, ringNormals, ringRadii, 1);
    }
 /*   std::vector<std::string> strs;
    if (it == 0)
    {
      std::cout << "drawing text" << std::endl;
      debug->drawText(ringCentres, strs);
    }*/
  }
  // finally we should filter these trunks to avoid overlapping ones
  std::cout << "num trunks: " << trunks.size() << std::endl;
  
  std::vector<Trunk> newTrunks;
  for (auto &trunk1: trunks)
  {
    bool isBest = true;
    for (auto &trunk2: trunks)
    {
      if (&trunk1==&trunk2)
        continue;
      double dist = vector2d(trunk1.centre - trunk2.centre).squaredNorm();
      if ((dist < sqr(trunk1.radius + trunk2.radius)) && trunk2.score > trunk1.score)
      {
        isBest = false;
        break;
      }
    }
    if (isBest)
      newTrunks.push_back(trunk1);
  }
  trunks = newTrunks;
  std::cout << "num trunks after pruning: " << trunks.size() << std::endl;
  
  if (verbose)
  {
    DebugDraw::instance()->drawTrunks(trunks);
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }
}
#else
#include "raygrid.h"
#include "raycuboid.h"

static const double trunk_height_to_width = 3.0; // height extent relative to real radius of trunk
static const double boundary_radius_scale = 2.0; // how much farther out is the expected boundary compared to real trunk radius? Larger requires more space to declare it a trunk

void getOverlap(const Grid<Eigen::Vector3d> &grid, const Trunk &trunk, std::vector<Eigen::Vector3d> &points)
{
  Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
  Eigen::Vector3d base = trunk.centre - 0.5*trunk.length*lean;
  Eigen::Vector3d top = trunk.centre + 0.5*trunk.length*lean;
  double outer_radius = trunk.radius * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, 0);
  Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
//  std::cout << "radius: " << trunk.radius << ", min: " << cuboid.min_bound_.transpose() << ", max: " << cuboid.max_bound_.transpose() << std::endl;
//  std::cout << "box min: " << grid.box_min.transpose() << std::endl;
//  std::cout << "mins: " << mins.transpose() << ", maxs: " << maxs.transpose() << std::endl;
  // just in case...
  mins = maxVector(mins, Eigen::Vector3i(0,0,0));
  Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1,1,1);
  maxs = minVector(maxs, min_dims);

  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0]<=maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1]<=maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2]<=maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        for (auto &pos: cell.data)
        {
          double h = pos[2] - trunk.centre[2];
          if (std::abs(h) > trunk.length*0.5)
          {
            continue;
          }
          Eigen::Vector3d centre = trunk.centre + lean * h;
          double dist2 = (centre - pos).squaredNorm();
          if (dist2 <= outer_radius*outer_radius)
          {
            points.push_back(pos);
          }
        }
      }
    }
  }
}


void getOverlap(const Grid<Trunk> &grid, const Trunk &trunk, std::vector<Trunk*> &overlapping_trunks)
{
  Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
  const double length_scale = 2.0; // make trunks longer when searching for overlapping trunks
  Eigen::Vector3d base = trunk.centre - 0.5*trunk.length*length_scale*lean;
  Eigen::Vector3d top = trunk.centre + 0.5*trunk.length*length_scale*lean;
  double outer_radius = trunk.radius * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, 0);
  Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();

  mins = maxVector(mins, Eigen::Vector3i(0,0,0));
  Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1,1,1);
  maxs = minVector(maxs, min_dims);

  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0]<=maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1]<=maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2]<=maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        for (auto &other_trunk: cell.data)
        {
          if (&other_trunk == &trunk)
            continue;
          // cylinder-cylinder overlap
          double lowest = std::max(trunk.centre[2]-0.5*trunk.length*length_scale, other_trunk.centre[2]-0.5*other_trunk.length*length_scale);
          double highest = std::min(trunk.centre[2]+0.5*trunk.length*length_scale, other_trunk.centre[2]+0.5*other_trunk.length*length_scale);
          if (highest <= lowest)
            continue;
          #if 1 // better
          Eigen::Vector3d oth = other_trunk.centre + vector3d(other_trunk.lean, 1)*(trunk.centre[2]-other_trunk.centre[2]) - trunk.centre;
          Eigen::Vector3d oth_dir = vector3d(other_trunk.lean, 1) - vector3d(trunk.lean, 1);
          oth_dir[2] = 0;
          double d = std::max(lowest-trunk.centre[2], std::min(oth.dot(oth_dir)/oth_dir.dot(oth_dir), highest-trunk.centre[2]));
          Eigen::Vector3d intersection = oth + oth_dir*d;
          double dist = intersection.norm();
          #else
          double mid = (lowest + highest)/2.0;
          Eigen::Vector3d p1 = trunk.centre + vector3d(trunk.lean, 1)*(mid - trunk.centre[2]);
          Eigen::Vector3d p2 = other_trunk.centre + vector3d(other_trunk.lean, 1)*(mid - other_trunk.centre[2]);
          Eigen::Vector3d p = (p2 - p1);
          double dist = p.norm();
          #endif
          if (dist > 1.25*(trunk.radius + other_trunk.radius))
            continue;
          overlapping_trunks.push_back((Trunk *)&other_trunk);
        }
      }
    }
  }
}

Wood::Wood(const Cloud &cloud, double midRadius, double, bool verbose)
{
  double spacing = cloud.estimatePointSpacing();
  
  // Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
  // trunk thickness affects how strictly it adheres to a cylinder.
  // lower trunk_thickness (stricter) will require a lower minimum_score to find the same number of trees.
  // at the point where minimum score is 0, it is invariant to the number of points
  const double minimum_score = 0.2/sqr(spacing);
  const double trunk_thickness = 0.025; 
  
  if (verbose)
  {
    std::cout << "estimated point spacig: " << spacing << ", minimum score: " << minimum_score << std::endl;
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }

  Eigen::Vector3d min_bound = cloud.calcMinBound();
  Eigen::Vector3d max_bound = cloud.calcMaxBound();
  std::cout << "cloud from: " << min_bound.transpose() << " to: " << max_bound.transpose() << std::endl;

  // 1. voxel grid of points (an acceleration structure)
  const double voxel_width = midRadius * 2.0;
  Grid<Eigen::Vector3d> grid(min_bound, max_bound, voxel_width);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3d pos = cloud.ends[i];
    grid.insert(grid.index(pos), pos);
  }
  const int min_num_points = 10;

  // 2. initialise one trunk candidate for each occupied voxel
  std::vector<Trunk> trunks;
  grid.walkCells([&](const Grid<Eigen::Vector3d> &grid, const Grid<Eigen::Vector3d>::Cell &cell)
  {
    if (cell.data.size() < min_num_points) // not enough data to attempt
      return;
    Trunk trunk;
    trunk.centre = (cell.index.cast<double>() + Eigen::Vector3d(0.5,0.5,0.5)) * voxel_width + grid.box_min;
    trunk.radius = midRadius;
    trunk.length = 2.0 * midRadius * trunk_height_to_width;
    trunk.score = trunk.weight = 0;
    trunk.lean.setZero();
    trunks.push_back(trunk);
  });

  // 3. iterate every candidate several times
  const int num_iterations = 5;
  std::vector<double> final_scores;
  std::vector<Eigen::Vector3d> final_points;
  for (int it = 0; it<num_iterations; it++)
  {
    double avScore = 0;
    double num = 0;
    for (int trunk_id = 0; trunk_id < (int)trunks.size(); trunk_id++)
    {
      auto &trunk = trunks[trunk_id];
      // get overlapping points to this trunk
      std::vector<Eigen::Vector3d> points;
      getOverlap(grid, trunk, points);
      if (points.size() < min_num_points) // not enough data to use
      {
        trunks[trunk_id] = trunks.back(); // so remove the trunk
        trunks.pop_back();
        trunk_id--;
        continue;
      }

      // weight these points
      Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
      std::vector<double> weights(points.size());

      Accumulator sum;
      trunk.score = 0;
      std::vector<double> scores(points.size());
      std::vector<Eigen::Vector3d> ps(points.size());
      struct Acc
      {
        Acc(){ x2 = y2 = xy = xz = yz = 0; }
        double x2, y2, xy, xz, yz;
      };
      Acc plane;

      Eigen::Vector3d mean_p(0,0,0);
      for (size_t i = 0; i<points.size(); i++)
      {
        double h = points[i][2] - trunk.centre[2];
        Eigen::Vector2d offset = vector2d(points[i] - (trunk.centre + lean*h));
#define PARABOLOID_METHOD
#if defined PARABOLOID_METHOD
        Eigen::Vector2d xy = offset/trunk.radius;
        double l2 = xy.squaredNorm();
        Eigen::Vector3d point(xy[0], xy[1], 0.5*l2); // a paraboloid that has gradient 1 at 1
        ps[i] = point;
        mean_p += point;
#endif


        double dist = offset.norm();
        double w = 1.0 - dist/(trunk.radius * boundary_radius_scale); // lateral fade off
        weights[i] = w;
        // remove radius. If radius_removal_factor=0 then half-sided trees will have estimated trunk centred on that edge
        //                If radius_removal_factor=1 then v thin trunks may accidentally get a radius and it won't shrink down
        const double radius_removal_factor = 0.5;
        offset -= offset * radius_removal_factor * trunk.radius / offset.norm(); 



        // lean, shift and change radius
        sum.x += h*w;
        sum.y += offset*w;
        sum.xy += h*offset*w;
        sum.x2 += h*h*w;
        sum.radius += dist;
        sum.radius2 += std::abs(h)*w;
        sum.weight += w;      

        // hard coding for now. Representing the expected error from circular in metres for real trees
        double score_centre = 1.0 - trunk.radius/trunk_thickness;
        double score_radius = 1.0;
        double score_2radius = 1.0 - trunk.radius/trunk_thickness;
        double score_3radius = 1.0 - trunk.radius/trunk_thickness;
        double weight = 0.0;
        if (dist < trunk.radius)
          weight = score_centre + (score_radius - score_centre)*dist/trunk.radius;
        else if (dist < 2.0*trunk.radius)
          weight = score_radius + (score_2radius - score_radius)*(dist - trunk.radius)/trunk.radius;
        else
          weight = score_2radius + (score_3radius - score_2radius)*(dist - 2.0*trunk.radius)/trunk.radius;

        if (it == num_iterations-1)
        {
          scores[i] = weight;
        }
        trunk.score += weight;
      }
      mean_p /= (double)points.size();
      for (auto &p: ps)
      {
        Eigen::Vector3d q = p - mean_p;
        plane.x2 += q[0]*q[0];
        plane.y2 += q[1]*q[1];
        plane.xy += q[0]*q[1];        
        plane.xz += q[0]*q[2];        
        plane.yz += q[1]*q[2];        
      }
      double n = sum.weight;
      trunk.score /= 2.0 * kPi * trunk.radius * trunk.length; // normalize

      avScore += trunk.score;
      num++;

      if (it == num_iterations-1 && trunk.score > minimum_score)
      {
        final_scores.insert(final_scores.end(), scores.begin(), scores.end());
        final_points.insert(final_points.end(), points.begin(), points.end());        
      }
      if (it == num_iterations-1 && trunk.score < minimum_score) // then remove the trunk
      {
        trunks[trunk_id] = trunks.back(); 
        trunks.pop_back();
        trunk_id--;
        continue;        
      }

      // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
      Eigen::Vector2d sXY = sum.xy - sum.x*sum.y/n;
      double sXX = sum.x2 - sum.x*sum.x/n;

      trunk.lean += sXY / sXX;
      double l = trunk.lean.norm();
      const double max_lean = 0.5;//25;
      if (l > max_lean)
        trunk.lean *= max_lean/l;          
      
#if defined PARABOLOID_METHOD
      double A = (plane.xz*plane.y2 - plane.yz*plane.xy) / (plane.x2*plane.y2 - plane.xy*plane.xy);
      double B = (plane.yz - A * plane.xy) / plane.y2;

      Eigen::Vector2d shift(A,B);
      double height = mean_p[2] + (shift[0]-mean_p[0])*A + (shift[1]-mean_p[1])*B;
      double paraboloid_height = 0.5*(A*A + B*B);
      double new_radius = std::sqrt(2.0*(height - paraboloid_height)) * trunk.radius;

      new_radius = std::min(new_radius, trunk.radius * 2.0); // don't grow by more than a factor of 2 per iteration.
      double shift2 = shift.squaredNorm();
      if (shift2 > 1.0) // don't shift more than one radius each iteration
        shift /= std::sqrt(shift2);

      trunk.centre += vector3d(shift * trunk.radius);   

      double radius_scale = new_radius / trunk.radius;
 //     double radius_scale = (sum.radius/(double)points.size()) / trunk.radius;
#else
      double radius_scale = (sum.radius/(double)points.size()) / trunk.radius;
      trunk.centre += vector3d(sum.y/n);
#endif
      trunk.centre[2] += sum.x / n;
      double length_scale = (sum.radius2/n) / (trunk.length * 0.25);
//      double scale = (radius_scale + length_scale)/2.0; // average, since we're not affecting the ratio of length to radius
//      trunk.radius *= scale;
//      trunk.length *= scale;
      trunk.radius *= radius_scale;
      trunk.length *= length_scale;
      if (trunk.radius > 0.5*trunk.length || trunk.length < midRadius) // not enough data to use
      {
        trunks[trunk_id] = trunks.back(); // so remove the trunk
        trunks.pop_back();
        trunk_id--;
        continue;
      }
    }
    if (verbose)
    {
      std::cout << "average score: " << avScore / num << std::endl;
      DebugDraw::instance()->drawTrunks(trunks);
    }  
  }
  if (verbose)
  {
 //   DebugDraw::instance()->drawCloud(final_points, final_scores, 1);
    DebugDraw::instance()->drawTrunks(trunks);
  } 

  // now I need to connect all the trunks into tree shapes
  Grid<Trunk> trunk_grid(min_bound, max_bound, voxel_width);
  for (auto &trunk: trunks)
  {
    trunk.next_down = NULL;
    trunk_grid.insert(trunk_grid.index(trunk.centre), trunk);
  }
  for (size_t i = 0; i<trunks.size(); i++)
  {
    std::vector<Trunk*> overlaps;
    getOverlap(trunk_grid, trunks[i], overlaps);
    for (auto &overlap: overlaps)
    {
      // no, not quite, we need an insertion sort type thing...
      if (trunks[i].centre[2] > overlap->centre[2])
      {
        if (trunks[i].next_down == NULL || trunks[i].next_down->centre[2] < overlap->centre[2])
          trunks[i].next_down = overlap;
      }
      else
      {
        if (overlap->next_down == NULL || overlap->next_down->centre[2] < trunks[i].centre[2])
          overlap->next_down = &trunks[i];
      }
    }
  }
  for (auto &trunk: trunks) // TODO: we can also walk this and generate a tree of branches
  {
    Trunk *tr = &trunk;
    while (tr->next_down != NULL)
      tr = tr->next_down;
    if (tr->radius != -1)
      trunk_bases.push_back(*tr);
    tr->radius = -1;
  }
  double mean_radius = 0;
  double total_volume = 0;
  const double volume_gain = 35.0; // this is the factor for the artificial trees in raycreate
  for (auto &t: trunk_bases)
  {
    mean_radius += t.radius;
    total_volume += volume_gain * t.radius * t.radius * t.radius;
  }
  mean_radius /= (double)trunk_bases.size();
  std::cout << "number of trees found: " << trunk_bases.size() << " with mean diameter: " << 2.0*mean_radius << std::endl;
  std::cout << "approximate volume of wood: " << total_volume << " cubic metres" << std::endl;
  const double wood_density = 500; // for Red gum, see: https://www.engineeringtoolbox.com/wood-density-d_40.html
  std::cout << "estimated mass of trees (at 500kg/m^3): " << total_volume * wood_density / 1000.0 << " tonnes" << std::endl;
  if (verbose)
  {
    DebugDraw::instance()->drawTrunks(trunk_bases);
  } 
}

bool Wood::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "Tree base location list: x, y, z, radius" << std::endl;
  for (auto &trunk: trunk_bases)
  {
    Eigen::Vector3d base = trunk.centre - vector3d(trunk.lean, 1)*trunk.length*0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << trunk.radius << std::endl;
  }
  return true;
}


#endif
} // namespace ray
