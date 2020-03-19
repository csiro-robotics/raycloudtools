// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raydraw.h"
#include <nabo/nabo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << endl;
  cout << "usage:" << endl;
  cout << "rayalign raycloudA raycloudB." << endl;
  exit(error);
}

struct Covariance
{
  Vector3d pos;
  Vector3d vectors[3];
  Vector3d values;
  double score;

  double getWidth() // trying to get the widest horizontal width of the ellipsoid
  {
    // do it as a weighted average for each vector pair.
    double weights[3];
    double widths[3];
    for (int i = 0; i<3; i++)
    {
      int j = (i+1)%3;
      widths[i] = max(values[i], values[j]);
      weights[i] = 1.0/(1e-4 + max(abs(vectors[i][2]), abs(vectors[j][2])));
    }

    return (widths[0]*weights[0] + widths[1]*weights[1] + widths[2]*weights[2])/(weights[0] + weights[1] + weights[2]);
  }
};

double crossCorrelate(const vector<double> &p1, const vector<double> &p2)
{
  // I'll do least squares for now, then later do reweighted least squares for better robustness
  return mean(p2) - mean(p1);
}


void getClosestVectors(const MatrixXd &pointsQ, const MatrixXd &pointsP, MatrixXi &indices, MatrixXd &dist2, int maxNeighbours, double maxDistance)
{
  cout << "get closest vectors. Vector length: " << pointsQ.rows() << ", number of points: " << pointsQ.cols() << endl;
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  if (maxNeighbours == 1)
    nns = Nabo::NNSearchD::createBruteForce(pointsP, pointsQ.rows());
  else
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, pointsQ.rows());
  indices.resize(maxNeighbours, pointsQ.cols());
  dist2.resize(maxNeighbours, pointsQ.cols());
  nns->knn(pointsQ, indices, dist2, maxNeighbours, 0.01, Nabo::NNSearchD::SORT_RESULTS, maxDistance);
  delete nns;
  cout << "get closest vectors complete" << endl;
}


int main(int argc, char *argv[])
{
  ros::init(argc, argv, "rayalign");
  DebugDraw draw;
  if (argc != 3)
    usage();

  string fileA = argv[1];
  string fileB = argv[2];
  Cloud cloudA, cloudB;
  cloudA.load(fileA);
  cloudB.load(fileB);

  Cloud *cloud[2] = {&cloudA, &cloudB};
  if (cloudA.ends.size() > cloudB.ends.size())
    swap(cloud[0], cloud[1]);
  bool swapped = cloud[0] == &cloudB;

  // This is a tricky algorithm, it should go like so:
  // 1. decimate
  // 2. get 20 nearest neighbours
  // 3. iterate each point, calculate its ellipsoid then add it to one of three sorted lists:
  //    a. direction list: these should be maximally vertical, thin and the nearest 20 centroid should be close to the nearest 15 (i.e. invariant to angle)
  //    b. ground/ceiling list: these should be maximally horizontal
  //    c. keypoint list: these should maximise the difference between 20 centroid and 15 centroid, i.e. corner-like
  // 4. crop these lists each to some percentage of the number of points (perhaps 10%)
  // 5. get the closest direction-list ellipsoids in the larger cloud to those in the smaller cloud, based on ellipsoid shape
  // 6. iterate through the list for candidate pairs, for each pair:
  //    a. transform the larger cloud onto the smaller cloud based on this position and direction
  //    b. use multiple vertical e^ix on the nearest ground/ceiling ellipsoids, to work out height offset, translate the cloud
  //    c. use multiple horizontal (tangent to direction ellipsoid) e^ix against key points to work out horizontal offset
  //    d. assess the amount of overlap between ellipsoids for this candidate
  // 7. after all pairs tested, pick the best fitting and transform cloudA

  vector<Covariance> wallLists[2], floorLists[2], cornerLists[2];
  for (int c = 0; c<2; c++)
  {
    vector<Vector3d> &ends = cloud[c]->ends;
    draw.drawCloud(ends, 1.0, 0);
    int searchSize = 20;
    Nabo::NNSearchD *nns;
    Nabo::Parameters params("bucketSize", 8);
    MatrixXd pointsP(3, ends.size());
    for (int i = 0; i<(int)ends.size(); i++)
      pointsP.col(i) = ends[i];
    cout << "creating linear heap" << endl;
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);

    // Run the search
    MatrixXi indices;
    MatrixXd dists2;
    indices.resize(searchSize, ends.size());
    dists2.resize(searchSize, ends.size());
    cout << "knn" << endl;
    nns->knn(pointsP, indices, dists2, searchSize, 0.01, 0, 1.0);
    delete nns;

    cout << "iterating: " << ends.size() << " points" << endl;
    for (int i = 0; i<(int)ends.size(); i++)
    {
      Matrix3d scatter;
      scatter.setZero();
      Vector3d centroid(0,0,0);
      int numNeighbours = 0;
      for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
      {
        centroid += ends[indices(j, i)];
        numNeighbours++;
      }   
      if (numNeighbours == 0)
        continue;
      centroid /= (double)numNeighbours;
      Vector3d smallCentroid(0,0,0);
      int half = max(1, numNeighbours/2);
      for (int j = 0; j<half; j++)
        smallCentroid += ends[indices(j, i)];
      smallCentroid /= half;    
      double centroidDelta = (smallCentroid - centroid).norm();
      for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
      {
        Vector3d offset = ends[indices(j,i)] - centroid;
        scatter += offset * offset.transpose();
      }
      scatter /= (double)numNeighbours;

      SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
      ASSERT(eigenSolver.info() == Success); 

      Vector3d eigenValue = eigenSolver.eigenvalues();
      const double eps = 1e-5;
      if (eigenValue[0] <= eps || eigenValue[1] <= eps || eigenValue[2] <= eps)
        continue;
      Matrix3d eigenVector = eigenSolver.eigenvectors();

      Covariance covariance;
      covariance.pos = centroid;
      covariance.vectors[0] = eigenVector.col(0);
      covariance.vectors[1] = eigenVector.col(1);
      covariance.vectors[2] = eigenVector.col(2);
      covariance.values = Vector3d(sqrt(eigenValue[0]), sqrt(eigenValue[1]), sqrt(eigenValue[2]));

      // these are the heuristics for being a wall, floor or corner
      double wallScore = 1.0/(eps + abs(covariance.vectors[0][2] * covariance.values[0] * centroidDelta));
      double floorScore = sqr(covariance.vectors[0][2]);
      double cornerScore = 1.0/covariance.getWidth(); // * covariance.values[0]/covariance.values[1];

      // now enter them in each list
      covariance.score = wallScore;
      wallLists[c].push_back(covariance);
      covariance.score = floorScore;
      floorLists[c].push_back(covariance);
      covariance.score = cornerScore;
      cornerLists[c].push_back(covariance);
    }
    struct 
    {
      bool operator()(const Covariance &a, const Covariance &b) const { return a.score > b.score; } 
    } higherScore;
    sort(wallLists[c].begin(), wallLists[c].end(), higherScore);
    sort(floorLists[c].begin(), floorLists[c].end(), higherScore);
    sort(cornerLists[c].begin(), cornerLists[c].end(), higherScore);
    cout << "wallScore 0: " << wallLists[c][0].score << ", 1: " << wallLists[c][1].score << ", end: " << wallLists[c].back().score << endl;
    cout << "floorScore 0: " << floorLists[c][0].score << ", 1: " << floorLists[c][1].score << ", end: " << floorLists[c].back().score << endl;

    double listRatio = 0.05;
    int listSize = (double)ends.size() * listRatio; 
    cout << "list size: " << listSize << endl;
    wallLists[c].resize(listSize);
    floorLists[c].resize(listSize);
    cornerLists[c].resize(listSize);

    // now draw it...
    Vector3d colours[3] = {Vector3d(1,0,0), Vector3d(0,1,0), Vector3d(0,0,1)};
    for (int l = 0; l<3; l++)
    {
      vector<Covariance> &list = l==0 ? wallLists[c] : (l==1 ? floorLists[c] : cornerLists[c]);
      vector<Vector3d> centres(listSize), radii(listSize);
      vector<Matrix3d> matrices(listSize);
      for (int i = 0; i<listSize; i++)
      {
        centres[i] = list[i].pos;
        matrices[i].col(0) = list[i].vectors[0];
        matrices[i].col(1) = list[i].vectors[1];
        matrices[i].col(2) = list[i].vectors[2];
        radii[i] = list[i].values;
      }
      draw.drawEllipsoids(centres, matrices, radii, colours[l], l+3*c);
    }
  }

  // next get the closest wallLists in cloud b to cloud a, by shape...

  // then find the candidate pairs
  vector<Vector3i> wallPairs;
  {
    // how do we parameterise these ellipsoids?
    // by normal[2]? by values?
    // Answer: by only values... normal[2] could be close to zero for most buildings
    // well, we could make it the distance at the top of the ellipsoid... so in the same space... OK.

    vector<Covariance> &list1 = wallLists[0];
    vector<Covariance> &list2 = wallLists[1];

    // TODO: shall we add an extra factor here? how about mean abs height diff?
    MatrixXd pointsQ(4, list1.size());
    for (int i = 0; i<(int)list1.size(); i++)
    {
      double tiltDistance = abs(list1[i].vectors[0][2])*list1[i].values.norm();
      pointsQ.col(i) = Vector4d(list1[i].values[0], list1[i].values[1], list1[i].values[2], tiltDistance);
    }
    MatrixXd pointsP(4, list2.size());
    for (int i = 0; i<(int)list2.size(); i++)
    {
      double tiltDistance = abs(list2[i].vectors[0][2])*list2[i].values.norm();
      pointsP.col(i) = Vector4d(list2[i].values[0], list2[i].values[1], list2[i].values[2], tiltDistance);
    }
    MatrixXi indices;
    MatrixXd dist2;
    getClosestVectors(pointsQ, pointsP, indices, dist2, 1, 1.0);

    // now we need to iterate through to see matching pairs
    for (int i = 0; i<(int)list1.size(); i++)
      if (indices(0,i)>-1)
        wallPairs.push_back(Vector3i(i, indices(0,i), 10000.0*sqrt(dist2(0,i))));
    cout << "finished generating wall pairs" << endl;
    struct 
    {
      bool operator()(const Vector3i &a, const Vector3i &b) const { return a[2] < b[2]; } 
    } lessDifference;
    cout << "sorting" << endl;
    sort(wallPairs.begin(), wallPairs.end(), lessDifference);
    cout << "finished sorting" << endl;
  }

  // now find the closest floor and corners to the wall ellipsoids:
  MatrixXi indices[2][2];
  int neighbourSize = 5;
  for (int t = 0; t<2; t++) // ellipsoid type (floor or corner)
  {
    if (t==0)
      cout << "finding closest floors to wall points" << endl;
    else
      cout << "finding closest corners to wall points" << endl;
    for (int c = 0; c<2; c++) // cloud 
    {
      vector<Covariance> &list1 = wallLists[c];
      vector<Covariance> &list2 = t==0 ? floorLists[c] : cornerLists[c];
      MatrixXd pointsQ(3, list1.size());
      for (int i = 0; i<(int)list1.size(); i++)
        pointsQ.col(i) = list1[i].pos;
      MatrixXd pointsP(3, list2.size());
      for (int i = 0; i<(int)list2.size(); i++)
        pointsP.col(i) = list2[i].pos;

      MatrixXd dist2;
      getClosestVectors(pointsQ, pointsP, indices[t][c], dist2, neighbourSize, 10.0);
    }    
  }

  // next, iterate through wallPairs
  double bestProximity = 0.0;
  Pose bestTransform;
  cout << "iterating through wall pairs, there are " << wallPairs.size() << endl;
  for (auto &pair: wallPairs)
  {
    vector<Vector3d> pairStarts, pairEnds;
    Pose poses[2];
    for (int i = 0; i<2; i++)
    {
      Vector3d &n = wallLists[i][pair[i]].vectors[0];
      poses[i].rotation = Quaterniond(AngleAxisd(atan2(n[0], n[1]), Vector3d(0,0,1)));
      poses[i].position = wallLists[i][pair[i]].pos;
      cout << "pose: " << poses[i] << endl;
    }
    pairStarts.push_back(wallLists[0][pair[0]].pos);
    pairEnds.push_back(wallLists[1][pair[1]].pos);

    // first, floors:
    {
      vector<double> points[2];
      for (int c = 0; c<2; c++)
      {
        MatrixXi &nearestFloors = indices[0][c];
        int i = pair[c];
        for (int j = 0; j<neighbourSize && nearestFloors(j,i)>-1; j++)
        {
          int id = nearestFloors(j,i);  
          points[c].push_back(floorLists[c][id].pos[2]);
          pairStarts.push_back(wallLists[c][pair[c]].pos);
          pairEnds.push_back(floorLists[c][id].pos);
        } 
      }
      double heightOffset = crossCorrelate(points[0], points[1]);
      cout << "height offset: " << heightOffset << endl;
      poses[1].position[2] = poses[0].position[2] + heightOffset;
    }

    // then corners:
    {
      vector<double> points[2];
      Vector3d sides[2];
      for (int c = 0; c<2; c++)
      {
        MatrixXi &nearestCorners = indices[1][c];
        int i = pair[c];
        sides[c] = wallLists[c][pair[c]].vectors[0].cross(Vector3d(0,0,1));
        for (int j = 0; j<neighbourSize && nearestCorners(j,i)>-1; j++)
        {
          int id = nearestCorners(j,i);  
          points[c].push_back((cornerLists[c][id].pos - wallLists[c][pair[c]].pos).dot(sides[c]));
          pairStarts.push_back(wallLists[c][pair[c]].pos);
          pairEnds.push_back(cornerLists[c][id].pos);
        } 
      }
      double sideOffset = crossCorrelate(points[0], points[1]);
      cout << "side offset: " << sideOffset << endl;
      poses[1].position += sideOffset*sides[1];
    }
    draw.drawLines(pairStarts, pairEnds);
    Pose transform = poses[1] * ~poses[0];
    cout << "now finding closest after applying the transform: " << transform << endl;
    // now rotate all the corners by the transform, and get a metric of how close they are, using closest points again!
    int list1Size = cornerLists[0].size();
    MatrixXd pointsQ(3, list1Size);
    for (int i = 0; i<list1Size; i++)
      pointsQ.col(i) = transform * cornerLists[0][i].pos; // transform the smaller cloud (faster)
    int list2Size = cornerLists[1].size();
    MatrixXd pointsP(3, list2Size);
    for (int i = 0; i<list2Size; i++)
      pointsP.col(i) = cornerLists[1][i].pos; 
    MatrixXi indices;
    MatrixXd dists2;
    getClosestVectors(pointsQ, pointsP, indices, dists2, 1, 1.0);
    double proximity = 0.0;
    for (int i = 0; i<list1Size; i++)
      if (indices(0,i)>-1)
        proximity += 1.0 / (0.1 + dists2(0,i));
    if (proximity > bestProximity)
    {
      bestProximity = proximity;
      bestTransform = transform;
    }
  }
  cout << "best proximity: " << bestProximity << " with transformation: " << bestTransform << endl;
  // finally, apply the transformation to the correct ray cloud and save it out. 
  if (swapped)
    bestTransform = ~bestTransform;
  cloudA.transform(bestTransform, 0.0);

  string fileStub = fileA;
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);
  cloudA.save(fileStub + "_aligned.ply");  

  return true;
}
