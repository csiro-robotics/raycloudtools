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
};

double crossCorrelate(const vector<double> &p1, const vector<double> &p2)
{
  return 0.0;
}


void getClosestVectors(const MatrixXd &pointsQ, const MatrixXd &pointsP, MatrixXi &indices, MatrixXd &dist2, int maxNeighbours, double maxDistance)
{
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
}


int main(int argc, char *argv[])
{
  ros::init(argc, argv, "rayalign");
  DebugDraw draw;
  if (argc != 3)
    usage();

  string fileA = argv[1];
  Cloud cloud[2];
  cloud[0].load(fileA);

  string fileB = argv[2];
  cloud[1].load(fileB);

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
    vector<Vector3d> &ends = cloud[c].ends;
    draw.drawCloud(ends, 1.0, 0);
    int searchSize = 20;
    Nabo::NNSearchD *nns;
    Nabo::Parameters params("bucketSize", 8);
    MatrixXd pointsP(3, ends.size());
    for (unsigned int i = 0; i<ends.size(); i++)
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
    for (unsigned int i = 0; i<(int)ends.size(); i++)
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
      double cornerScore = centroidDelta * covariance.values[0]/covariance.values[1];

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
    vector<Covariance> *lists[3] = {&wallLists[c], &floorLists[c], &cornerLists[c]};
    Vector3d colours[3] = {Vector3d(1,0,0), Vector3d(0,1,0), Vector3d(0,0,1)};
    for (int l = 0; l<3; l++)
    {
      vector<Covariance> &list = l==0 ? wallLists[c] : (l==1 ? floorLists[c] : cornerLists[c]);
      vector<Vector3d> centres(listSize), radii(listSize);
      vector<Matrix3d> matrices(listSize);
      vector<Vector3d> p2s(listSize);
      for (int i = 0; i<listSize; i++)
      {
        centres[i] = list[i].pos;
        matrices[i].col(0) = list[i].vectors[0];
        matrices[i].col(1) = list[i].vectors[1];
        matrices[i].col(2) = list[i].vectors[2];
        radii[i] = list[i].values;
        p2s[i] = list[i].pos + list[i].vectors[0];
      }
      draw.drawEllipsoids(centres, matrices, radii, colours[l], l);
  //    draw.drawLines(centres, p2s);
    }
    cin.get();
  }

  int a = cloud[0].ends.size() < cloud[1].ends.size() ? 0 : 1;
  int b = 1-a;
  // next get the closest wallLists in cloud b to cloud a, by shape...

  // then find the candidate pairs
  vector<Vector3i> wallPairs;
  {
    // how do we parameterise these ellipsoids?
    // by normal[2]? by values?
    // Answer: by only values... normal[2] could be close to zero for most buildings
    // well, we could make it the distance at the top of the ellipsoid... so in the same space... OK.

    vector<Covariance> &list1 = wallLists[a];
    vector<Covariance> &list2 = wallLists[b];

    int searchSize = 1;
    MatrixXd pointsQ(4, list1.size());
    for (unsigned int i = 0; i<list1.size(); i++)
    {
      double tiltDistance = abs(list1[i].vectors[0][2])*list1[i].values.norm();
      pointsQ.col(i) = Vector4d(list1[i].values[0], list1[i].values[1], list1[i].values[2], tiltDistance);
    }
    MatrixXd pointsP(4, list2.size());
    for (unsigned int i = 0; i<list2.size(); i++)
    {
      double tiltDistance = abs(list2[i].vectors[0][2])*list2[i].values.norm();
      pointsP.col(i) = Vector4d(list2[i].values[0], list2[i].values[1], list2[i].values[2], tiltDistance);
    }
    MatrixXi indices;
    MatrixXd dist2;
    getClosestVectors(pointsQ, pointsP, indices, dist2, 1, 1.0);

    // now we need to iterate through to see matching pairs
    for (unsigned int i = 0; i<list1.size(); i++)
    {
      vector<Vector3d> centres(2), radii(2);
      vector<Matrix3d> mats(2);

      centres[0] = list1[i].pos;
      mats[0].col(0) = list1[i].vectors[0];
      mats[0].col(1) = list1[i].vectors[1];
      mats[0].col(2) = list1[i].vectors[2];
      radii[0] = list1[i].values;

      if (indices(0,i)>-1)
      {
        // TODO: check that it is actually sorting the list
        int id = indices(0,i);

        centres[1] = list2[id].pos;
        mats[1].col(0) = list2[id].vectors[0];
        mats[1].col(1) = list2[id].vectors[1];
        mats[1].col(2) = list2[id].vectors[2];
        radii[1] = list2[id].values;
      
        draw.drawEllipsoids(centres, mats, radii, Vector3d(1,1,1), 0);

        // here's where we now use the floor and corners to guide us...
        // first we need to get the closest floor ellipsoids to each wall ellipsoid.
        wallPairs.push_back(Vector3i(i, id, 10000.0*sqrt(dist2(0,i))));
      }
  //    cin.get();
    }
    struct 
    {
      bool operator()(const Vector3i &a, const Vector3i &b) const { return a[3] < b[3]; } 
    } lessDifference;
    sort(wallPairs.begin(), wallPairs.end(), lessDifference);
  }

  // now find the closest floor and corners to the wall ellipsoids:
  MatrixXi indices[2][2];
  int neighbourSize = 5;
  for (int t = 0; t<2; t++) // ellipsoid type (floor or corner)
  {
    for (int l = 0; l<2; l++) // cloud 
    {
      int c = l==0 ? a : b;
      vector<Covariance> &list1 = wallLists[c];
      vector<Covariance> &list2 = t==0 ? floorLists[c] : cornerLists[c];
      MatrixXd pointsQ(3, list1.size());
      for (unsigned int i = 0; i<list1.size(); i++)
        pointsQ.col(i) = list1[i].pos;
      MatrixXd pointsP(3, list2.size());
      for (unsigned int i = 0; i<list2.size(); i++)
        pointsP.col(i) = list2[i].pos;

      MatrixXd dist2;
      getClosestVectors(pointsQ, pointsP, indices[t][c], dist2, neighbourSize, 10.0);
    }    
  }

  // next, iterate through wallPairs
  double bestProximity = 0.0;
  Vector3d bestTranslation(0,0,0);
  double bestRotation = 0.0;

  for (auto &pair: wallPairs)
  {
    Vector3d translation = wallLists[a][pair[0]].pos - wallLists[b][pair[0]].pos;
    Vector3d norm1 = wallLists[a][pair[0]].vectors[0];
    Vector3d norm2 = wallLists[b][pair[1]].vectors[0];
    norm1[2] = 0.0;
    norm2[2] = 0.0;
    double rotation = atan2(norm1.dot(norm2), norm1.cross(norm2)[2]);
    cout << "translation: " << translation.transpose() << ", rotation: " << rotation << endl;

    // first, floors:
    {
      vector<double> points[2];
      for (int l = 0; l<2; l++)
      {
        int c = l==0 ? a : b;
        MatrixXi &nearestFloors = indices[0][l];
        int i = pair[l];
        for (int j = 0; j<neighbourSize && nearestFloors(j,i)>-1; j++)
        {
          int id = nearestFloors(j,i);  
          points[c].push_back(floorLists[c][id].pos[2]);
        } 
      }
      double heightOffset = crossCorrelate(points[0], points[1]);
      cout << "height offset: " << heightOffset << endl;
      translation[2] += heightOffset;
    }

    // then corners:
    {
      vector<double> points[2];
      Vector3d sides[2];
      for (int l = 0; l<2; l++)
      {
        int c = l==0 ? a : b;
        MatrixXi &nearestCorners = indices[0][l];
        int i = pair[l];
        sides[l] = wallLists[c][pair[l]].vectors[0].cross(Vector3d(0,0,1));
        for (int j = 0; j<neighbourSize && nearestCorners(j,i)>-1; j++)
        {
          int id = nearestCorners(j,i);  
          points[c].push_back(cornerLists[c][id].pos.dot(sides[l]));
        } 
      }
      double sideOffset = crossCorrelate(points[0], points[1]);
      cout << "side offset: " << sideOffset << endl;
      translation += sideOffset*sides[0];
    }

    // now rotate all the corners by the transform, and get a metric of how close they are, using closest points again!


    int list1Size = cornerLists[a].size();
    MatrixXd pointsQ(3, list1Size);
    for (unsigned int i = 0; i<list1Size; i++)
      pointsQ.col(i) = cornerLists[a][i].pos;
    int list2Size = cornerLists[b].size();
    MatrixXd pointsP(3, list2Size);
    for (unsigned int i = 0; i<list2Size; i++)
      pointsP.col(i) = cornerLists[b][i].pos; // transform!
    MatrixXi indices;
    MatrixXd dists2;
    getClosestVectors(pointsQ, pointsP, indices, dists2, 1, 1.0);
    double proximity = 0.0;
    for (int i = 0; i<list1Size; i++)
      if (indices(0,i)>-1)
        proximity += 1.0 / (0.1 + dists2(0,i));
    if (proximity > bestProximity)
    {
      bestTranslation = translation;
      bestRotation = rotation;
    }
  }

  // finally, apply the transformation to the correct ray cloud and save it out. 
/*  Cloud newCloud = clouds[0];
  Pose pose(bestTranslation, Quaterniond(AngleAxisd(bestRotation, Vector3d(0,0,1))));
  newCloud.transform(pose, 0.0);

  string fileStub = files[0];
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);
  newCloud.save(fileStub + "_aligned.ply");  
*/
  return true;
}
