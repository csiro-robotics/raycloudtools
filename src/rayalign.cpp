// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "rayalignment.h"
#include "raydraw.h"
#include <nabo/nabo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <complex>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << endl;
  cout << "usage:" << endl;
  cout << "rayalign raycloudA raycloudB" << endl;
  cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << endl;
  exit(error);
}

void getMeanAndNormal(const vector<Vector3d> &points, const vector<int> &ids, Vector3d &centroid, Vector3d &normal, Vector3d &width, Matrix3d *mat=NULL, Vector3d *vec = NULL)
{
  Vector3d total(0,0,0);
  for (auto &id: ids)
    total += points[id];
  centroid = total / (double)ids.size();

  Matrix3d scatter;
  scatter.setZero();
  for (auto &id: ids)
  {
    Vector3d offset = points[id] - centroid;
    scatter += offset * offset.transpose();
  }
  scatter / (double)ids.size();

  SelfAdjointEigenSolver<Matrix3d> eigenSolver(scatter.transpose());
  ASSERT(eigenSolver.info() == Success); 
  normal = eigenSolver.eigenvectors().col(0);
  width = eigenSolver.eigenvalues();
  if (mat)
    *mat = eigenSolver.eigenvectors();
  if (vec)
  {
    *vec = eigenSolver.eigenvalues();
    *vec = maxVector(*vec, Vector3d(1e-5,1e-5,1e-5));
  }
}

int main(int argc, char *argv[])
{
  // TODO: This method works when there is more than 30% overlap. 
  // For less, we can additionally repeat this procedure for small scale cubic sections (perhaps 20% of map width)
  // then use the fft power spectrum (low res) as a descriptor into the knn to find matching locations, 
  // then pick the best match.
  #if defined(USE_ROS)
    ros::init(argc, argv, "rayalign");
  #endif
    DebugDraw draw;

  if (argc != 3 && argc != 4)
    usage();
  bool verbose = false;
  if (argc == 4)
  {
    if (string(argv[3]) != "--verbose" && string(argv[3]) != "-v")
      usage();
    verbose = true;
  }

  string fileA = argv[1];
  string fileB = argv[2];

  AlignTranslationYaw aligner;
  aligner.clouds[0].load(fileA);
  aligner.clouds[1].load(fileB);
  
  aligner.alignCloud0ToCloud1(0.5, verbose);
  string fileStub = fileA;
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);
  if (verbose)
    aligner.clouds[0].save(fileStub + "_coarse_aligned.ply");  

  // Next we need a fine grained alignment. 
  // Method: 
  // grid the points and pick one point from each grid cell
  // then (using the 7 neighbours), find the points within a radius of just these.... (or should I use NABO?)
  // then, get ellipsoid for each point.
  // finally, run it-reweighted-least-squares to fit to the bend model

  vector<Matrix3d> matrices[2];
  vector<Vector3d> radii[2];

  vector<Vector3d> centroids[2];
  vector<Vector3d> normals[2];
  vector<Vector3d> widths[2];
  for (int c = 0; c<2; c++)
  {
    // 1. decimate quite fine
    vector<int> decimated = voxelSubsample(aligner.clouds[c].ends, 0.1); 
    vector<Vector3d> decimatedPoints;
    decimatedPoints.reserve(decimated.size());
    vector<Vector3d> decimatedStarts;
    decimatedStarts.reserve(decimated.size());
    for (int i = 0; i<(int)decimated.size(); i++)
    {
      if (aligner.clouds[c].rayBounded(decimated[i]))
      {
        decimatedPoints.push_back(aligner.clouds[c].ends[decimated[i]]);
        decimatedStarts.push_back(aligner.clouds[c].starts[decimated[i]]);
      }
    }

    // 2. find the coarser random candidate points. We just want a fairly even spread but not the voxel centres
    vector<int> candidates = voxelSubsample(decimatedPoints, 1.0);
    vector<Vector3d> candidatePoints(candidates.size());
    vector<Vector3d> candidateStarts(candidates.size());
    for (int i = 0; i<(int)candidates.size(); i++)
    {
      candidatePoints[i] = decimatedPoints[candidates[i]];   
      candidateStarts[i] = decimatedStarts[candidates[i]];
    }

    int searchSize = 20;
    int qSize = candidates.size();
    int pSize = decimatedPoints.size();

    Nabo::NNSearchD *nns;
    MatrixXd pointsQ(3, qSize);
    for (int i = 0; i<qSize; i++)
      pointsQ.col(i) = candidatePoints[i];
    MatrixXd pointsP(3, pSize);
    for (int i = 0; i<pSize; i++)
      pointsP.col(i) = decimatedPoints[i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 3);

    // Run the search
    MatrixXi indices;
    MatrixXd dists2;
    indices.resize(searchSize, qSize);
    dists2.resize(searchSize, qSize);
    nns->knn(pointsQ, indices, dists2, searchSize, 0.01, 0, 1.0);
    delete nns;

    centroids[c].reserve(qSize);
    normals[c].reserve(qSize);
    widths[c].reserve(qSize);
    matrices[c].reserve(qSize);
    radii[c].reserve(qSize);
    vector<int> ids;
    ids.reserve(searchSize);
    const int minPointsPerEllipsoid = 5;
    for (int i = 0; i<qSize; i++)
    {
      ids.clear();
      for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
        ids.push_back(indices(j,i));
      if (ids.size() < minPointsPerEllipsoid) // not dense enough
        continue; 

      Vector3d normal;
      Vector3d centroid;
      Vector3d width;
      getMeanAndNormal(decimatedPoints, ids, centroid, normal, width);
      if ((centroid - candidateStarts[i]).dot(normal) > 0.0)
        normal = -normal;     
      // now repeat but removing back facing points. This deals better with double walls, which are quite common
      for (int j = (int)ids.size()-1; j>=0; j--)
      {
        int id = ids[j];
        if ((decimatedPoints[id] - decimatedStarts[id]).dot(normal) > 0.0)
        {
          ids[j] = ids.back();
          ids.pop_back();
        }
      }
      if (ids.size() < minPointsPerEllipsoid) // not dense enough
        continue; 
      Matrix3d mat;
      Vector3d radius;
      getMeanAndNormal(decimatedPoints, ids, centroid, normal, width, &mat, &radius);
      if (width[0]/width[1] > sqr(0.5)) // not planar enough
        continue;
      if ((centroid - candidateStarts[i]).dot(normal) > 0.0)
        normal = -normal;     
      matrices[c].push_back(mat);
      radii[c].push_back(radius);
      centroids[c].push_back(centroid);
      normals[c].push_back(normal);       
      widths[c].push_back(width);  
    }
    Vector3d cols[2] = {Vector3d(1,0,0), Vector3d(0,1,0)};
    draw.drawEllipsoids(centroids[c], matrices[c], radii[c], cols[c], c);
  }
  
  // Now we have the ellipsoids, we need to do a nearest neighbour on this set of data, trying to match orientation as well as location
  {
    vector<Vector3d> lineStarts;
    vector<Vector3d> lineEnds;
    int searchSize = 1;
    int qSize = centroids[0].size();
    int pSize = centroids[1].size();

    Nabo::NNSearchD *nns;
    MatrixXd pointsQ(6, qSize);
    for (int i = 0; i<qSize; i++)
      pointsQ.col(i) << centroids[0][i], normals[0][i];
    MatrixXd pointsP(6, pSize);
    for (int i = 0; i<pSize; i++)
      pointsP.col(i) << centroids[1][i], normals[1][i];
    nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 6);

    // Run the search
    double maxRange = 1.0;
    MatrixXi indices;
    MatrixXd dists2;
    indices.resize(searchSize, qSize);
    dists2.resize(searchSize, qSize);
    nns->knn(pointsQ, indices, dists2, searchSize, 0.01, 0, maxRange);
    delete nns;

    struct Match
    {
      Vector3d pos[2];
      Vector3d normal;
      double weight;
    };
    vector<Match> matches;
    vector<int> ids;
    ids.reserve(searchSize);
    for (int i = 0; i<qSize; i++)
    {
      // shall we pick only two-way matches? Not for now...
      for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
      {
        lineStarts.push_back(centroids[0][i]);
        lineEnds.push_back(centroids[1][indices(j,i)]);
        Match match;
        match.pos[0] = centroids[0][i];
        match.pos[1] = centroids[1][indices(j,i)];
        match.normal = (normals[0][i] + normals[1][indices(j,i)]).normalized();
        match.weight = 1.0/(0.03 + sqrt(widths[0][i][0] + widths[1][i][0]));
        matches.push_back(match);
      }
    }
    draw.drawLines(lineStarts, lineEnds);
    // OK, what else can I do?
    // 1. two-way matches
    // 2. weight based on thinness  // doesn't seem to help!
    // 3. reweight based on orientation difference too
    // 4. match cylinders differently
    // 5. match spheres differently
    // 6. or at least ignore non-planars // doesn't seem to help in this case, but left in at 0.5 eccentricity

    // with these matches we now run the iterative reweighted least squares..
    Pose pose = Pose::identity();
    int maxIterations = 5;
    for (int it = 0; it<maxIterations; it++)
    {
      double d = 10.0*(double)it/(double)maxIterations;
      // stage 1: only Euclidean. 
      const int stateSize = 6;
      Matrix<double, stateSize, stateSize> AtA;
      AtA.setZero();
      Matrix<double, stateSize, 1> AtB;
      AtB.setZero();
      double squareError = 0.0;
      for (auto &match: matches)
      {
        double error = (match.pos[1] - match.pos[0]).dot(match.normal);
        double weight = /*match.weight **/ pow(max(1.0 - sqr(error/maxRange), 0.0), d*d);
        squareError += sqr(error);
        Matrix<double, stateSize, 1> At;

        for (int i = 0; i<3; i++) // change in error with change in raycloud translation
          At[i] = match.normal[i];
        for (int i = 0; i<3; i++) // change in error with change in raycloud orientation
        {
          Vector3d axis(0,0,0);
          axis[i] = 1.0;
          At[3+i] = -(match.pos[0].cross(axis)).dot(match.normal);
        }
        AtA += At*weight*At.transpose();
        AtB += At*weight*error;
      }
      Matrix<double, stateSize, 1> x = AtA.ldlt().solve(AtB);
      cout << "rmse: " << sqrt(squareError/(double)matches.size()) << endl;
      cout << "least squares shift: " << x[0] << ", " << x[1] << ", " << x[2] << ", rotation: " << x[3] << ", " << x[4] << ", " << x[5] << endl;
      Vector3d rot(x[3], x[4], x[5]);
      double angle = rot.norm();
      rot.normalize();
      Quaterniond halfRot(AngleAxisd(angle/2.0, rot));
      Pose shift(Vector3d(x[0], x[1], x[2]), Quaterniond(AngleAxisd(angle, rot)));
      pose = shift * pose;
      for (auto &match: matches)
      {
        match.pos[0] = shift*match.pos[0];
        match.normal = halfRot * match.normal;
      }
      cout << "pose: " << pose << endl;
    }
    aligner.clouds[0].transform(pose, 0.0);
  }

  aligner.clouds[0].save(fileStub + "_aligned.ply");  

  return true;
}
