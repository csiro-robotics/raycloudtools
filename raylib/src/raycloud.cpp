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

void Cloud::getSurfels(int searchSize, vector<Vector3d> *centroids, vector<Vector3d> *normals, vector<Vector3d> *dimensions, vector<Matrix3d> *mats, MatrixXi *neighbourIndices)
{
  // simplest scheme... find 3 nearest neighbours and do cross product
  if (centroids)
    centroids->resize(ends.size());
  if (normals)
    normals->resize(ends.size());
  if (dimensions)
    dimensions->resize(ends.size());
  if (mats)
    mats->resize(ends.size());
  Nabo::NNSearchD *nns;
  vector<int> rayIDs;
  rayIDs.reserve(ends.size());
  for (unsigned int i = 0; i<ends.size(); i++)
    if (rayBounded(i))
      rayIDs.push_back(i);
  MatrixXd pointsP(3, rayIDs.size());
  for (unsigned int i = 0; i<rayIDs.size(); i++)
    pointsP.col(i) = ends[rayIDs[i]];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);
  
  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(searchSize, rayIDs.size());
  dists2.resize(searchSize, rayIDs.size());
  nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0); // TODO: needs to sort here
  delete nns;

  if (neighbourIndices)
    neighbourIndices->resize(searchSize, ends.size());
  for (int i = 0; i<(int)rayIDs.size(); i++)
  {
    int I = rayIDs[i];
    if (neighbourIndices)
    {
      int j;
      for (j = 0; j<searchSize && indices(j,i)>-1; j++)
        (*neighbourIndices)(j,I) = rayIDs[indices(j,i)];    
      if (j<searchSize)
        (*neighbourIndices)(j,I) = -1;
    }

    Vector3d centroid = ends[I];
    int num;
    for (num = 0; num<searchSize && indices(num,i)>-1; num++)
      centroid += ends[rayIDs[indices(num,i)]];
    centroid /= (double)(num+1);
    if (centroids)
      (*centroids)[I] = centroid;

    Matrix3d scatter = (ends[I] - centroid) * (ends[I] - centroid).transpose();
    for (int j = 0; j<num; j++)
    {
      Vector3d offset = ends[rayIDs[indices(j,i)]] - centroid;
      scatter += offset * offset.transpose();
    }
    scatter /= (double)(num + 1);

    SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
    ASSERT(eigenSolver.info() == Success);
    if (normals)
    { 
      Vector3d normal = eigenSolver.eigenvectors().col(0);
      if ((ends[I] - starts[I]).dot(normal) > 0.0)
        normal = -normal;
      (*normals)[I] = normal;
    }
    if (dimensions)
    {
      Vector3d eigenvals = maxVector(Vector3d(1e-10,1e-10,1e-10), eigenSolver.eigenvalues());
      (*dimensions)[I] = Vector3d(sqrt(eigenvals[0]), sqrt(eigenvals[1]), sqrt(eigenvals[2]));
    }
    if (mats)
      (*mats)[I] = eigenSolver.eigenvectors();
  }
}

// starts are required to get the normal the right way around
vector<Vector3d> Cloud::generateNormals(int searchSize)
{
  vector<Vector3d> normals;
  getSurfels(searchSize, NULL, &normals, NULL, NULL, NULL);
  return normals;
}

// TODO: this could call getSurfels too!
void Cloud::generateEllipsoids(vector<Ellipsoid> &ellipsoids)
{
  cout << "generating " << ends.size() << " ellipsoids" << endl;
  ellipsoids.resize(ends.size());
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
    ellipsoids[i].transient = false;
    if (!rayBounded(i))
    {
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
  }
}

void fillGrid(Grid<int> &grid, const vector<Vector3d> &starts, const vector<Vector3d> &ends)
{
  cout << "filling grid with " << ends.size() << " rays" << endl;
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

void Cloud::markIntersectedEllipsoids(Grid<int> &grid, vector<bool> &transients, vector<Ellipsoid> &ellipsoids, const string &mergeType, double numRays, bool selfTransient)
{
  cout << "mark intersected ellipsoids, numRays: " << numRays << ", mergeType: " << mergeType << endl;
  DebugDraw draw;
  draw.drawCloud(ends, 1.0, 0);
  int type = mergeType == "oldest" ? 0 : (mergeType == "newest" ? 1 : (mergeType == "min" ? 2 : 3));
  vector<bool> rayTested;
  rayTested.resize(ends.size(), false);
  int cnt = 0;
  for (auto &ellipsoid: ellipsoids)
  {
    if ((cnt++)%20000 == 0)
      cout << cnt << "/" << ellipsoids.size() << endl;

//    if (numRays>0 && (ellipsoid.pos - Vector3d(2,-0.5,0.5)).norm() < 0.2)
//      cout << "point" << endl;
    if (ellipsoid.transient) // a previous ellipsoid could have removed the ray that represents this ellipsoid
      continue; 
    if (ellipsoid.extents == Vector3d::Zero()) // unbounded rays cannot be a transient object
      continue;
    // get all the rays that overlap this ellipsoid
    Vector3d bMin = (ellipsoid.pos - ellipsoid.extents - grid.boxMin)/grid.voxelWidth;
    Vector3d bMax = (ellipsoid.pos + ellipsoid.extents - grid.boxMin)/grid.voxelWidth;
    if (bMax[0] < 0.0 || bMax[1] < 0.0 || bMax[2] < 0.0)
      continue;
    if (bMin[0] > (double)grid.dims[0]-1.0 || bMin[1] > (double)grid.dims[1]-1.0 || bMin[2] > (double)grid.dims[2]-1.0)
      continue;
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
      rayTested[rayID] = false;

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
    double numBefore = 0, numAfter = 0;
    ellipsoid.numRays = hits + passThroughIDs.size();
    if (numRays == 0 || selfTransient)
      ellipsoid.opacity = (double)hits / ((double)hits + (double)passThroughIDs.size());
    if (ellipsoid.numRays == 0 || ellipsoid.opacity == 0 || numRays == 0) 
      continue;
    if (selfTransient) 
    {
      ellipsoid.numGone = passThroughIDs.size();
      // now get some density stats...
      double misses = 0;
      for (auto &i: passThroughIDs)
      {
        if (times[i] > lastIntersectionTime)
          numAfter++;
        else if (times[i] < firstIntersectionTime)
          numBefore++;
        else
          misses++;
      }
      double h = hits + 1e-8 - 1.0; // subtracting 1 gives an unbiased opacity estimate
      ellipsoid.opacity = h/(h + misses);
      ellipsoid.numGone = numBefore + numAfter;
    }
    else // compare to other cloud
    {
      if (passThroughIDs.size() > 0)
      {
        if (times[passThroughIDs[0]] > ellipsoid.time)
          numAfter = (double)passThroughIDs.size();
        else 
          numBefore = (double)passThroughIDs.size();
      }
    }
    
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
      double d = 0.0;
      for (int j = 0; j<(int)passThroughIDs.size(); j++)
      {
        d += ellipsoid.opacity;
        if (d>=1.0)
          d--;
        else
          continue;
        int i = passThroughIDs[j];
        if (!selfTransient || times[i] < firstIntersectionTime || times[i] > lastIntersectionTime)
        {
          // remove ray i
          transients[i] = true;
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
  
  vector<Ellipsoid> ellipsoids;
  generateEllipsoids(ellipsoids);

  Grid<int> grid(boxMin, boxMax, voxelWidth);
  fillGrid(grid, starts, ends);

  vector<bool> transients;
  transients.resize(ends.size(), false);
  // now walk every ray through the grid and mark if transient
  markIntersectedEllipsoids(grid, transients, ellipsoids, mergeType, numRays, true);

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
    if (ellipsoids[i].transient || transients[i])
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
  vector<Grid<int>> grids(clouds.size());
  for (int c = 0; c<(int)clouds.size(); c++)
  {
    Vector3d boxMin(1e10,1e10,1e10);
    Vector3d boxMax(-1e10,-1e10,-1e10);
    for (auto &end: clouds[c].ends)
    {
      boxMin = minVector(boxMin, end);
      boxMax = maxVector(boxMax, end);
    }
    grids[c].init(boxMin, boxMax, voxelWidth);
    fillGrid(grids[c], clouds[c].starts, clouds[c].ends);
  }

  vector<vector<bool> > transients(clouds.size());
  for (int c = 0; c<(int)clouds.size(); c++)
    transients[c].resize(clouds[c].ends.size(), false);
  // now for each cloud, look for other clouds that penetrate it
  for (int c = 0; c<(int)clouds.size(); c++)
  {
    vector<Ellipsoid> ellipsoids;
    clouds[c].generateEllipsoids(ellipsoids);
    // just set opacity
    clouds[c].markIntersectedEllipsoids(grids[c], transients[c], ellipsoids, mergeType, 0, false);
    
    for (int d = 0; d<(int)clouds.size(); d++)
    {
      if (d==c)
        continue;
      // use ellipsoid opacity to set transient flag true on transients
      clouds[d].markIntersectedEllipsoids(grids[d], transients[d], ellipsoids, mergeType, numRays, false);
    }

    for (int i = 0; i<(int)clouds[c].ends.size(); i++)
      if (ellipsoids[i].transient)
        transients[c][i] = true; // HACK, we need a better way to signal this!
  }
  for (int c = 0; c<(int)clouds.size(); c++)
  {
    auto &cloud = clouds[c];
    int t = 0;
    int f = 0;
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      if (transients[c][i])
      {
        t++;
        differences.starts.push_back(cloud.starts[i]);
        differences.ends.push_back(cloud.ends[i]);
        differences.times.push_back(cloud.times[i]);
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

void Cloud::split(Cloud &cloud1, Cloud &cloud2, function<bool(int i)> fptr)
{
  for (int i = 0; i<(int)ends.size(); i++)
  {
    Cloud &cloud = fptr(i) ? cloud2 : cloud1;
    cloud.starts.push_back(starts[i]);
    cloud.ends.push_back(ends[i]);
    cloud.times.push_back(times[i]);
    cloud.colours.push_back(colours[i]);
  }
}
