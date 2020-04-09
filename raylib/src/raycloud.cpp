// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
#include "raydraw.h"
#include <nabo/nabo.h>
#include <set>

using namespace std;
using namespace Eigen;
using namespace RAY;

void Cloud::save(const std::string &fileName)
{
  string name = fileName;
  if (name.substr(name.length()-4) != ".ply")
    name += ".ply";
  writePly(name, starts, ends, times, colours);
}

bool Cloud::load(const std::string &fileName)
{
  // look first for the raycloud PLY
  if (fileName.substr(fileName.size()-4)==".ply")
    return loadPLY(fileName);
  if(ifstream((fileName + ".ply").c_str(), ios::in))
    return loadPLY(fileName + ".ply");

  // otherwise, look for a .laz and _traj.txt file by that name
  if (ifstream((fileName + ".laz").c_str(), ios::in) && ifstream((fileName + "_traj.txt").c_str(), ios::in))
    return loadLazTraj(fileName + ".laz", fileName + "_traj.txt");

  return false;
}  

bool Cloud::load(const std::string &pointCloud, const std::string &trajFile)
{
  string nameEnd = pointCloud.substr(pointCloud.size()-4);
  if (nameEnd == ".ply")
    readPly(pointCloud, starts, ends, times, colours);
  else if (nameEnd == ".laz" || nameEnd == ".las")
    readLas(pointCloud, ends, times, colours, 1);   
  else
  {
    cout << "Error converting unknown type: " << pointCloud << endl;
    return false;
  }

  Trajectory trajectory;
  trajectory.load(trajFile);
  calculateStarts(trajectory);
  return true;  
}

bool Cloud::loadPLY(const string &file)
{
  return readPly(file, starts, ends, times, colours);
}

bool Cloud::loadLazTraj(const string &lazFile, const string &trajFile)
{
  bool success = readLas(lazFile, ends, times, colours, 1);
  if (!success)
    return false;
  Trajectory trajectory;
  trajectory.load(trajFile);
  calculateStarts(trajectory);
  return true;
}

void Cloud::calculateStarts(const Trajectory &trajectory)
{
  // Aha!, problem in calculating starts when times are not ordered.
  if (trajectory.nodes.size()>0)
  {
    int n = 1;
    starts.resize(ends.size());
    for (int i = 0; i<(int)ends.size(); i++)
    {
      while ((times[i] > trajectory.nodes[n].time) && n<(int)trajectory.nodes.size()-1)
        n++;
      double blend = (times[i] - trajectory.nodes[n-1].time)/(trajectory.nodes[n].time-trajectory.nodes[n-1].time);
      starts[i] = trajectory.nodes[n-1].pose.position + (trajectory.nodes[n].pose.position - trajectory.nodes[n-1].pose.position)*clamped(blend, 0.0, 1.0);
    }
  }
  else
    cout << "can only recalculate when a trajectory is available" << endl;
}

Vector3d Cloud::calcMinBound()
{
  Vector3d minV(1e10, 1e10, 1e10);
  for (int i = 0; i<(int)ends.size(); i++)
  {
    if (rayBounded(i))
      minV = minVector(minV, ends[i]);
  }
  return minV;
}

Vector3d Cloud::calcMaxBound()
{
  Vector3d maxV(-1e10, -1e10, -1e10);
  for (int i = 0; i<(int)ends.size(); i++)
  {
    if (rayBounded(i))
      maxV = maxVector(maxV, ends[i]);
  }
  return maxV;
}

void Cloud::transform(const Pose &pose, double timeDelta)
{
  for (int i = 0; i<(int)starts.size(); i++)
  {
    starts[i] = pose * starts[i];
    ends[i] = pose * ends[i];
    times[i] += timeDelta;
  }  
}

void Cloud::removeUnboundedRays()
{
  vector<int> valids;
  for (int i = 0; i<(int)ends.size(); i++)
    if (rayBounded(i))
      valids.push_back(i);
  for (int i = 0; i<(int)valids.size(); i++)
  {
    starts[i] = starts[valids[i]];
    ends[i] = ends[valids[i]];
    times[i] = times[valids[i]];
    colours[i] = colours[valids[i]];
  }
  starts.resize(valids.size());
  ends.resize(valids.size());
  times.resize(valids.size());
  colours.resize(valids.size());
}

void Cloud::decimate(double voxelWidth)
{
  vector<int> subsample = voxelSubsample(ends, voxelWidth);
  for (int i = 0; i<(int)subsample.size(); i++)
  {
    int id = subsample[i];
    starts[i] = starts[id];
    ends[i] = ends[id];
    colours[i] = colours[id];
    times[i] = times[id];
  }
  starts.resize(subsample.size());
  ends.resize(subsample.size());
  colours.resize(subsample.size());
  times.resize(subsample.size());
}

// starts are required to get the normal the right way around
vector<Vector3d> Cloud::generateNormals(int searchSize)
{
  // simplest scheme... find 3 nearest neighbours and do cross product
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  MatrixXd pointsP(3, ends.size());
  for (unsigned int i = 0; i<ends.size(); i++)
    pointsP.col(i) = ends[i];//col.transpose();
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);//, 0, params);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(searchSize, ends.size());
  dists2.resize(searchSize, ends.size());
  nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0);
  delete nns;

  vector<Vector3d> normals;
  for (unsigned int i = 0; i<ends.size(); i++)
  {
    Matrix3d scatter;
    scatter.setZero();
    Vector3d average(0,0,0);
    double numNeighbours = 0;
    for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
    {
      average += ends[indices(j, i)];
      numNeighbours++;
    }
    average /= numNeighbours;
    for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
    {
      Vector3d offset = ends[indices(j,i)] - average;
      scatter += offset * offset.transpose();
    }

    SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
    ASSERT(eigenSolver.info() == Success); 
    Matrix3d eigenVector = eigenSolver.eigenvectors();
    Vector3d normal = eigenVector.col(0);
    if ((ends[i] - starts[i]).dot(normal) > 0.0)
      normal = -normal;
    normals.push_back(normal);
  }
  return normals;
}

void Cloud::generateEllipsoids(vector<Ellipsoid> &ellipsoids)
{
  int searchSize = 16;
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  MatrixXd pointsP(3, ends.size());
  for (unsigned int i = 0; i<ends.size(); i++)
    pointsP.col(i) = ends[i];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(searchSize, ends.size());
  dists2.resize(searchSize, ends.size());
  nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0);
  delete nns;

  for (unsigned int i = 0; i<ends.size(); i++)
  {
    if (!rayBounded(i))
    {
      ellipsoids[i].transient = false;
      ellipsoids[i].extents.setZero();
      continue;
    }
    Matrix3d scatter;
    scatter.setZero();
    Vector3d centroid(0,0,0);
    double numNeighbours = 1e-10;
    for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
    {
      int index = indices(j, i);
      if (rayBounded(index))
      {
        centroid += ends[index];
        numNeighbours++;
      }
    }
    centroid /= numNeighbours;
    for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
    {
      int index = indices(j, i);
      if (rayBounded(index))
      {
        Vector3d offset = ends[index] - centroid;
        scatter += offset * offset.transpose();
      }
    }
    scatter /= numNeighbours;

    SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
    ASSERT(eigenSolver.info() == Success); 

    Vector3d eigenValue = eigenSolver.eigenvalues();
    Matrix3d eigenVector = eigenSolver.eigenvectors();
    
    ellipsoids[i].pos = centroid;
    double scale = 1.7; // this scale roughly matches the dimensions of a uniformly dense ellipsoid
    eigenValue[0] = scale*sqrt(max(1e-10,eigenValue[0]));
    eigenValue[1] = scale*sqrt(max(1e-10,eigenValue[1]));
    eigenValue[2] = scale*sqrt(max(1e-10,eigenValue[2]));
    ellipsoids[i].vectors[0] = eigenVector.col(0)/eigenValue[0];
    ellipsoids[i].vectors[1] = eigenVector.col(1)/eigenValue[1];
    ellipsoids[i].vectors[2] = eigenVector.col(2)/eigenValue[2];
    ellipsoids[i].time = times[i];
    ellipsoids[i].setExtents(eigenVector, eigenValue);
    ellipsoids[i].setPlanarity(eigenValue);
    ellipsoids[i].transient = false;
  }
}

void fillGrid(Grid<int> &grid, const vector<Vector3d> &starts, const vector<Vector3d> &ends)
{
  // next populate the grid with these ellipsoid centres
  for (int i = 0; i<(int)ends.size(); i++)
  {
    if (!(i%20000))
      cout << i << "/" << ends.size() << endl;
    Vector3d dir = ends[i] - starts[i];
    Vector3d dirSign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (starts[i] - grid.boxMin)/grid.voxelWidth;
    Vector3d end = (ends[i] - grid.boxMin)/grid.voxelWidth;
    Vector3i startIndex(start[0], start[1], start[2]);
    Vector3i endIndex(end[0], end[1], end[2]);
    double lengthSqr = (endIndex - startIndex).squaredNorm();
    Vector3i index = startIndex;
    for (;;)
    {
      grid.insert(index[0], index[1], index[2], i);

      if (index == endIndex || (index - startIndex).squaredNorm()>lengthSqr)
        break;
      Vector3d mid = grid.boxMin + grid.voxelWidth*Vector3d(index[0]+0.5, index[1]+0.5, index[2]+0.5);
      Vector3d nextBoundary = mid + 0.5*grid.voxelWidth*dirSign;
      Vector3d delta = nextBoundary - starts[i];
      Vector3d d(delta[0]/dir[0], delta[1]/dir[1], delta[2]/dir[2]);
      if (d[0] < d[1] && d[0]<d[2])
        index[0] += dirSign[0];
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += dirSign[1];
      else
        index[2] += dirSign[2];
    }
  }

  grid.report();
}

void Cloud::markIntersectedEllipsoids(Grid<int> &grid, vector<Ellipsoid> &ellipsoids, const string &mergeType, double numRays)
{
  cout << "mark intersected ellipsoids, numRays: " << numRays << ", mergeType: " << mergeType << endl;
  DebugDraw draw;
  draw.drawCloud(ends, 1.0, 0);
  int type = mergeType == "oldest" ? 0 : (mergeType == "newest" ? 1 : (mergeType == "min" ? 2 : 3));
  vector<bool> rayTested(ellipsoids.size());
  fill(begin(rayTested), end(rayTested), false);
  int cnt = 0;
  for (auto &ellipsoid: ellipsoids)
  {
    if ((cnt++)%20000 == 0)
      cout << cnt << "/" << ellipsoids.size() << endl;
    if (ellipsoid.transient) // a previous ellipsoid could have removed the ray that represents this ellipsoid
      continue; 
    if (ellipsoid.extents == Vector3d::Zero()) // unbounded rays cannot be a transient object
      continue;
    ellipsoid.transient = false;
    // get all the rays that overlap this ellipsoid
    Vector3d bMin = (ellipsoid.pos - ellipsoid.extents - grid.boxMin)/grid.voxelWidth;
    Vector3d bMax = (ellipsoid.pos + ellipsoid.extents - grid.boxMin)/grid.voxelWidth;
    Vector3i bmin = maxVector(Vector3i(0,0,0), Vector3i(bMin[0], bMin[1], bMin[2]));
    Vector3i bmax = minVector(Vector3i(bMax[0], bMax[1], bMax[2]), Vector3i(grid.dims[0]-1, grid.dims[1]-1, grid.dims[2]-1));
    
    vector<int> rayIDs;
    for (int x = bmin[0]; x<=bmax[0]; x++)
    {
      for (int y = bmin[1]; y<=bmax[1]; y++)
      {
        for (int z = bmin[2]; z<=bmax[2]; z++)
        {
          auto &list = grid.cell(x,y,z).data;
          for (auto &i: list)
          { 
            if (rayTested[i])
              continue;
            rayTested[i] = true;
            rayIDs.push_back(i);
          }        
        }
      }
    }
    for (auto &rayID: rayIDs)
      rayTested[rayID] = 0;

    double firstIntersectionTime = 1e10;
    double lastIntersectionTime = -1e10;
    int hits = 0;
    vector<double> passThroughIDs;
    for (auto &i: rayIDs)
    {
      Vector3d dir = ends[i] - starts[i];
      const double passDistance = 0.05;
      double ratio = passDistance / dir.norm();
      // ray-ellipsoid intersection
      Vector3d toSphere = ellipsoid.pos - starts[i];
      Vector3d ray;
      ray[0] = dir.dot(ellipsoid.vectors[0]);
      ray[1] = dir.dot(ellipsoid.vectors[1]);
      ray[2] = dir.dot(ellipsoid.vectors[2]);
      double rayLength = ray.norm();
      ray /= rayLength;
      Vector3d to;
      to[0] = toSphere.dot(ellipsoid.vectors[0]);
      to[1] = toSphere.dot(ellipsoid.vectors[1]);
      to[2] = toSphere.dot(ellipsoid.vectors[2]);

      double d = to.dot(ray);
      double dist2 = (to - ray*d).squaredNorm();

      if (dist2 > 1.0) // misses the ellipsoid
        continue;
      double alongDist = sqrt(1.0 - dist2);
      if (rayLength < d - alongDist) // doesn't reach the ellipsoid
        continue; 
      bool passThrough = rayLength*(1.0-ratio) > d + alongDist; // last number requires rays to pass some way past the object
      if (passThrough)
        passThroughIDs.push_back(i);
      else
      {
        hits++;
        firstIntersectionTime = min(firstIntersectionTime, times[i]);
        lastIntersectionTime = max(lastIntersectionTime, times[i]);
      }
    }
    ellipsoid.numRays = hits + passThroughIDs.size();
    ellipsoid.opacity = 0.0;
    ellipsoid.numGone = passThroughIDs.size();
    if (hits == 0) // if nothing hits this ellipsoid, that's pretty weird, though technically it might just be possible
      continue;
    // now get some density stats...
    double misses = 0, numBefore = 0, numAfter = 0;
    for (auto &i: passThroughIDs)
    {
      double time = times[i];
      if (time > lastIntersectionTime)
        numAfter++;
      else if (time < firstIntersectionTime)
        numBefore++;
      else
        misses++;
    }
    double h = hits + 1e-8 - 1.0; // subtracting 1 gives an unbiased opacity estimate
    ellipsoid.opacity = h/(h + misses);
    ellipsoid.numGone = numBefore + numAfter;
    
    double sequenceLength = numRays/ellipsoid.opacity;
    int removeEllipsoid = false;
    if (type == 0 || type == 1) 
    {
      if (max(numBefore, numAfter) < sequenceLength) 
        continue;
      if (type == 0) // oldest
        removeEllipsoid = numBefore >= sequenceLength; // if false then remove numAfter rays if > seqLength
      else if (type == 1) // newest
        removeEllipsoid = numAfter >= sequenceLength; // if false then remove numBefore rays if > seqLength
    }
    else
    {
      // we use sum rather than max below, because it better picks out moving objects that may have some
      // pass through rays before and after the hit points.
      if (numBefore + numAfter < sequenceLength) 
        continue;
      removeEllipsoid = type == 2; // min is remove ellipsoid, max is remove ray
    }

    if (removeEllipsoid) 
      ellipsoid.transient = true;
    else // if we don't remove the ellipsoid then we should remove numBefore and numAfter rays if they're greater than sequence length
    {
      int c = misses;
      for (int j = 0; j<(int)passThroughIDs.size(); j++)
      {
        if (c > hits-1)
        {
          c -= hits-1;
          continue; // we're removing in proportion to translucency of the ellipsoid we are keeping! Remove line to treat ellipsoid as solid
        }
        c += misses;
        int i = passThroughIDs[j];
        double time = times[i];
        if (time < firstIntersectionTime || time > lastIntersectionTime)
        {
          // remove ray i
          rayTested[i] = true; // prevents it being used again
          ellipsoids[i].transient = true;
        }
      }
    }
  }
}

static double voxelWidth = 0.25;

void Cloud::findTransients(Cloud &transient, Cloud &fixed, const string &mergeType, double numRays, bool colourCloud)
{
  cout << "find transients" << endl;
  Vector3d boxMin(1e10,1e10,1e10);
  Vector3d boxMax(-1e10,-1e10,-1e10);
  for (auto &end: ends)
  {
    boxMin = minVector(boxMin, end);
    boxMax = maxVector(boxMax, end);
  }
  
  vector<Ellipsoid> ellipsoids(ends.size());
  generateEllipsoids(ellipsoids);

  cout << "created normals and spheres" << endl;

  Grid<int> grid(boxMin, boxMax, voxelWidth);
  fillGrid(grid, starts, ends);

  cout << "grid is populated" << endl;

  // now walk every ray through the grid and mark if transient
  markIntersectedEllipsoids(grid, ellipsoids, mergeType, numRays);

  cout << "generating two clouds from flags" << endl;

  // Lastly, generate the new ray clouds from this sphere information
  for (int i = 0; i<(int)ellipsoids.size(); i++)
  {
    RGBA col = colours[i];
    if (colourCloud)
    {
      col.red = (uint8_t)((1.0 - ellipsoids[i].planarity) * 255.0);
      col.blue = (uint8_t)(ellipsoids[i].opacity * 255.0);
      col.green = (uint8_t)((double)ellipsoids[i].numGone/((double)ellipsoids[i].numGone + 10.0) * 255.0);
    }
    if (ellipsoids[i].transient)
    {
      transient.starts.push_back(starts[i]);
      transient.ends.push_back(ends[i]);
      transient.times.push_back(times[i]);
      transient.colours.push_back(col);
    }
    else
    {
      fixed.starts.push_back(starts[i]);
      fixed.ends.push_back(ends[i]);
      fixed.times.push_back(times[i]);
      fixed.colours.push_back(col);
    }
  }
}

void Cloud::combine(std::vector<Cloud> &clouds, Cloud &differences, const string &mergeType, double numRays)
{
  Vector3d boxMin(1e10,1e10,1e10);
  Vector3d boxMax(-1e10,-1e10,-1e10);
  for (auto &cloud: clouds)
  {
    for (auto &end: cloud.ends)
    {
      boxMin = minVector(boxMin, end);
      boxMax = maxVector(boxMax, end);
    }
  }
  // now for each cloud, look for other clouds that penetrate it
  for (int c = 0; c<(int)clouds.size(); c++)
  {
    Cloud &cloud = clouds[c];
  
    // first generate all ellipsoids and grids containing them, for each cloud
    vector<Ellipsoid> ellipsoids;
    ellipsoids.resize(cloud.ends.size());
    cout << "generating ellipsoids" << endl;
    cloud.generateEllipsoids(ellipsoids);
    Grid<int> grid(boxMin, boxMax, voxelWidth);
    fillGrid(grid, cloud.starts, cloud.ends);
    
    for (int d = 0; d<(int)clouds.size(); d++)
    {
      if (d==c)
        continue;
      Cloud &other = clouds[d];
      cout << "marking intersected ellipsoids" << endl;
      other.markIntersectedEllipsoids(grid, ellipsoids, mergeType, numRays);
    }

    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (ellipsoids[i].transient)
        cloud.times[i] = -abs(cloud.times[i]); // HACK, we need a better way to signal this!
    }
  }
  for (auto &cloud: clouds)
  {
    int t = 0;
    int f = 0;
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (cloud.times[i] < 0.0)
      {
        t++;
        differences.starts.push_back(cloud.starts[i]);
        differences.ends.push_back(cloud.ends[i]);
        differences.times.push_back(-cloud.times[i]);
        differences.colours.push_back(cloud.colours[i]);
      }
      else
      {
        f++;
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        colours.push_back(cloud.colours[i]);
      }
    }      
    cout << t << " transients, " << f << " fixed rays." << endl;
  }
}
