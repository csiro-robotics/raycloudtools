#include "wood.h"
#include "debugDraw.h"

static DebugDraw *debug;
#define LEAN

Vector2d vector2d(const Vector3d &v)
{
  return Vector2d(v[0], v[1]);
}
Vector3d vector3d(const Vector2d &v)
{
  return Vector3d(v[0], v[1], 0.0);
}

Wood::Wood(const RayCloud &cloud, double midRadius, double heightRange, const vector<Vector2d> &markedPoints)
{
  debug = new DebugDraw;
  debug->drawCloud(cloud.ends, 0, 0.5);
 // debug->drawRayCloud(cloud.starts, cloud.ends, 1);
  
  vector<Vector3d> markedCentre(markedPoints.size());
  vector<Vector3d> markedNormals(markedPoints.size());
  vector<double> markedRadii(markedPoints.size());
  
  const int minPointsPerTrunk = 15;
  const double minScore = 1.1;
  double maxRadius = midRadius * 2.0;

  width = (1.5 + 1e-5) * maxRadius; 
  Vector2i minCell = Vector2i(floor(cloud.minBound[0] / width), floor(cloud.minBound[1] / width));
  Vector2i maxCell(floor(cloud.maxBound[0] / width), floor(cloud.maxBound[1] / width));
  Vector2d minBound = Vector2d(minCell[0], minCell[1]) * width; 
  size = Vector2i(maxCell[0] + 1 - minCell[0], maxCell[1] + 1 - minCell[1]);
  grid.resize(size[0] * size[1]);
  for (int x = 0; x < size[0]; x++)
    for (int y = 0; y < size[1]; y++)
      grid[x + size[0] * y].minBound = minBound + Vector2d(x, y) * width;

  // populate the cells
  for (int i = 0; i < (int)cloud.ends.size(); i++)
  {
    Vector2i index(floor((cloud.ends[i][0] - minBound[0]) / width), floor((cloud.ends[i][1] - minBound[1]) / width));
    Cell &cell = grid[index[0] + size[0] * index[1]];
    cell.rays.push_back(Ray(cloud.starts[i], cloud.ends[i]));
    cell.height += cloud.ends[i][2];
  }

  Trunk trunk;
  trunk.radius = midRadius;
  trunk.score = 0;
  for (auto &cell : grid)
  {
    cell.height /= (1e-10 + (double)cell.rays.size());
    if (cell.rays.size() >= minPointsPerTrunk)
    {
      trunk.centre = vector3d(cell.minBound + width * Vector2d(0.5,0.5));
      trunk.length = heightRange / 3.0; // the catchment starts at twice this, so catches 2/3 of the passed in cloud
      trunk.centre[2] = cell.height;
      trunk.lean = Vector2d(0,0);
      trunk.thickness = 0.0;
      trunks.push_back(trunk);
    }
  }
  
  for (int i = 0; i<(int)markedCentre.size(); i++)
  {
    markedCentre[i] = vector3d(markedPoints[i]);
    Vector2i index(floor((markedCentre[i][0] - minBound[0]) / width), floor((markedCentre[i][1] - minBound[1]) / width));
    index[0] = clamped(index[0], 0, size[0]-1);
    index[1] = clamped(index[1], 0, size[1]-1);
    markedCentre[i][2] = grid[index[0] + size[0] * index[1]].height;
    if (markedCentre[i][2] == 0.0 && i>0)
      markedCentre[i][2] = markedCentre[i-1][2];
    markedNormals[i] = Vector3d(0,0,1);
    markedRadii[i] = 0.2;
  }
 // debug->drawRings(markedCentre, markedNormals, markedRadii, 0);
  debug->drawText(components(trunks, centre), vector<string>());

  
  const int maxIterations = 10;
  for (int it = 0; it < maxIterations; it++)
  {
    debug->drawTrunks(trunks);
    vector<Vector3d> spheres;
    vector<Vector3d> ringCentres;
    vector<Vector3d> ringNormals;
    vector<double> ringRadii;
    double bestScore = 0.0;
    double radiusScale = 2.0 - 0.6*(double)it/(double)(maxIterations-1); // pull in the average, to get a better estimation after it converges
    double lengthScale = 2.0 - (double)it/(double)maxIterations;
    for (int tr = (int)trunks.size()-1; tr>=0; tr--)
    {
      Trunk &trunk = trunks[tr];
      // for each trunk candidate, find the (max) 4 cells it overlaps 
      double outerRadius = trunk.radius * radiusScale;
      Vector2d minPos = vector2d(trunk.centre) - outerRadius*Vector2d(1,1) - minBound;
      Vector2d maxPos = vector2d(trunk.centre) + outerRadius*Vector2d(1,1) - minBound;
      Vector2i indexMin(max(0.0, floor(minPos[0] / width)), max(0.0, floor(minPos[1] / width)));
      Vector2i indexMax(min((int)floor(maxPos[0] / width), size[0]-1), min((int)floor(maxPos[1] / width), size[1]-1));
      Vector3d average(0, 0, 0);
      trunk.weight = 0;
      vector<Vector3d> ps;
      const bool accountForInwardRays = it>0;
      vector<Ray> rays[5];

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
          minh = min(minh, cell.height);
          for (auto &rawRay : cell.rays)
          {
            Ray ray = rawRay;
#if defined(LEAN) // unlean time
            Vector3d offset = vector3d(trunk.lean) * (ray.pos[2] - trunk.centre[2]); 
            ray.pos -= offset;
            ray.start -= offset; 
#endif
            Vector2d xy = vector2d(ray.pos - trunk.centre)/trunk.radius;
            double l2 = xy.squaredNorm();
            if (l2 > sqr(outerRadius / trunk.radius))
              continue;
            if (ray.pos[2] > heights[0] && ray.pos[2] < heights[5])
            {
              // we're transforming around a Riemann sphere at trunk.centre with 90 degrees at trunk.radius
              Vector3d point(2.0*xy[0] / (1.0 + l2), 2.0*xy[1] / (1.0 + l2), (1.0 - l2) / (1.0 + l2));
              // TODO: test using SVD instead of PCA (eigen decomp)
              ps.push_back(point);
              spheres.push_back((point + Vector3d(0,0,0.5)) * trunk.radius + trunk.centre);
              
              average += point;
              trunk.weight++;
              int h = 0; 
              while (ray.pos[2] > heights[h+1])
                h++;
              if (h >= 5)
                cout << "huh?" << endl;
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
      
      Matrix3d scatter;
      scatter.setZero();
      for (auto &p : ps)
        scatter += (p-average) * (p-average).transpose();

      SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter); 
      ASSERT(eigenSolver.info() == Success); // use EigenSolver if it isn't always self-adjoint
      Vector3d eigenValues = eigenSolver.eigenvalues();
      Matrix3d eigenVectors = eigenSolver.eigenvectors();

      Vector3d normal = eigenVectors.col(0); // least eigenvaluestruct Ray
      if (normal.dot(average) < 0.0)
        normal = -normal;
      ringNormals.push_back(normal);
      Vector2d flatNormal = vector2d(normal);

      // now work out centre and radius from projected stats
      double phi = atan2(flatNormal.norm(), normal[2]);
      double h = average.dot(normal);
      ringCentres.push_back((normal*h + Vector3d(0,0,0.5)) * trunk.radius + trunk.centre);
      ringRadii.push_back(sqrt(1.0 - sqr(h)) * trunk.radius);
      double angle = acos(h);
      
      // make sure circle doesn't cross infinity (-z) line.
      double a = phi - angle;
      if (a < -pi)
        a += 2.0*pi;
      double b = phi + angle;
      if (b > pi)
        b -= 2.0*pi;
      
      double lMin = tan(min(a, b)/2.0) * trunk.radius;
      double lMax = tan(max(a, b)/2.0) * trunk.radius;
      double length = (lMax + lMin) / 2.0;
      flatNormal.normalize();
      double newRadius = (lMax - lMin) / 2.0;
      
      double intersectionArea;
      double radDiff = abs(newRadius - trunk.radius);
      // from https://www.xarg.org/2016/07/calculate-the-intersection-area-of-two-circles/
      if (abs(newRadius - trunk.radius) < 1e-5 && length < 1e-5)
        intersectionArea = pi * sqr((newRadius + trunk.radius)/2.0);
      else if (length < radDiff)
        intersectionArea = pi * sqr(min(newRadius, trunk.radius));
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
        
        if (intersectionArea > pi * min(R2, r2) + 1e-9)
          cout << "weird intersection area" << endl;
      }
      double scoreScale = (pi*sqr(trunk.radius) * pi*sqr(newRadius))/(1e-10 + sqr(intersectionArea)); // this grows with the square of the distance between the disks
      const double pruningLateness = 1.0; // if you want it to prune more early and aggressively, try reducing this to 0.5 
      scoreScale = pow(scoreScale, pruningLateness);
      if (it == maxIterations-1)
        scoreScale = 1.0; // no advantage once the iterations are over
      
      trunk.centre += vector3d(flatNormal*length);
      trunk.radius = newRadius;
      
      vector<Ray> rays2;
     
      Accumulator cumulative[6];
      
      // now adjust the length and position
      for (int i = 1; i<6; i++)
      {
        cumulative[i] = cumulative[i-1];
        for (auto &ray: rays[i-1])
        {
          Vector3d dif = ray.pos - trunk.centre;
          Vector2d vec = vector2d(dif);
          double d = vec.norm();
          double dR = d - trunk.radius;
          Vector2d pos = vec * dR/d;

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
      Vector2d bestAbsLean(0,0);
      Vector2d bestOffset(0,0);
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
          Vector2d sXY = sum.xy - sum.x*sum.y/N;
          double sXX = sum.x2 - sum.x*sum.x/N;

          Vector2d lean = sXY / sXX;
          
          Vector2d absLean = lean + trunk.lean;
          Vector2d oldLean = lean;
          double l = absLean.norm();
          double maxLean = 0.25;
          if (l > maxLean)
            absLean *= maxLean/l;          
          lean = absLean - trunk.lean;
          
          Vector2d offset = sum.y/N - lean*sum.x/N;
          
          double thickness = sYY - lean.squaredNorm() * sXX;
          
#define PARALLEL_RADIUS // This adjusts the radius each segment, to get a smaller thickness. Otherwise the radius is adjusted after the segment is found
#if defined PARALLEL_RADIUS
          double dRadius = sum.radius/N + sum.z.dot(offset)/N + sum.xz.dot(lean)/N;
          thickness -= sqr(dRadius)*N;
          double radius = trunk.radius + dRadius; // we also quickly update the radius here, since length and lean may have changed
#endif      
          if (thickness < 0.0)
            cout << "weird, negative variance- " << oldLean.transpose() << endl;

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
        Vector3d dif = ray.pos - trunk.centre;
        Vector2d vec = vector2d(dif) - trunk.lean * dif[2];
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
          Vector2d toStart = vector2d(ray.start - ray.pos);
          Vector2d toCentre = vector2d(trunk.centre - ray.pos);
          if (toStart.dot(toCentre) > 0.0)
          {
            double dist1 = toCentre.norm();
            double dot = clamped(toCentre.dot(toStart) / toStart.squaredNorm(), 0.0, 1.0);
            double dist2 = (-toCentre + toStart * dot).norm();
            double delta = sqr(trunk.radius - dist2) - sqr(trunk.radius - dist1); // adjust the sum of square errors
            thickness += max(0.0, 4.0 * delta);
          }
        }
      }
      trunk.thickness = sqrt(thickness / trunk.weight);
      // TODO: try replacing max(eigen...) with sqr(radius) * trunk.weight. Should preference partial trees
      trunk.score = sqrt((heightVariance / trunk.length) * max(eigenValues[1], eigenValues[2]) / thickness); // max would do better with trees seen from one side, but more false positives
      // remove poorly performing trunks
      if (trunk.weight < minPointsPerTrunk || trunk.score*scoreScale < minScore)
      {
        trunks[tr] = trunks.back();
        trunks.pop_back();
        continue;
      }
      bestScore = max(bestScore, trunk.score);
    }
    cout << "best score: " << bestScore << endl;
    debug->drawCloud(spheres, 1, 1.0);
    debug->drawRings(ringCentres, ringNormals, ringRadii, 1);
    vector<string> strs;
  /*  if (it == 0)
    {
      cout << "drawing text" << endl;
      debug->drawText(ringCentres, strs);
    }*/
  }
  // finally we should filter these trunks to avoid overlapping ones
  cout << "num trunks: " << trunks.size() << endl;
  
  vector<Trunk> newTrunks;
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
  cout << "num trunks after pruning: " << trunks.size() << endl;
  
  debug->drawTrunks(trunks);
}
