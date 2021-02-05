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
Eigen::Vector3d vector3d(const Eigen::Vector2d &v)
{
  return Eigen::Vector3d(v[0], v[1], 0.0);
}

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
    DebugDraw::instance()->drawTrunks(trunks);
}
} // namespace ray
