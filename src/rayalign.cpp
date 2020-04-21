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
  cout << "                             --rigid    - rigid alignment only" << endl;
  cout << "                             --verbose  - outputs FFT images and the coarse alignment cloud" << endl;
  cout << "                             --local    - fine alignment only, assumes clouds are already approximately aligned" << endl;
  exit(error);
}

void getSurfel(const vector<Vector3d> &points, const vector<int> &ids, Vector3d &centroid, Vector3d &width, Matrix3d &mat)
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
  width = maxVector(eigenSolver.eigenvalues(), Vector3d(1e-5,1e-5,1e-5));
  width = Vector3d(sqrt(width[0]), sqrt(width[1]), sqrt(width[2]));
  mat = eigenSolver.eigenvectors();
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

  if (argc < 3 || argc > 6)
    usage();
  bool verbose = false;
  bool rigidOnly = false;
  bool localOnly = false;
  for (int a = 3; a<argc; a++)
  {
    if (string(argv[a]) == "--verbose" || string(argv[a]) == "-v")
      verbose = true;
    else if (string(argv[a]) == "--rigid" || string(argv[a]) == "-r")
      rigidOnly = true;
    else if (string(argv[a]) == "--local" || string(argv[a]) == "-l")
      localOnly = true;
    else
      usage();
  }

  string fileA = argv[1];
  string fileB = argv[2];
  string fileStub = fileA;
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);

  AlignTranslationYaw aligner;
  aligner.clouds[0].load(fileA);
  aligner.clouds[1].load(fileB);

  if (!localOnly)
  {
    aligner.alignCloud0ToCloud1(0.5, verbose);
    if (verbose)
      aligner.clouds[0].save(fileStub + "_coarse_aligned.ply");  
  }

  // Next the fine grained alignment. 
  // Method: 
  // 1. for each cloud: decimate the cloud to make it even, but still quite detailed, e.g. one point per cubic 10cm
  // 2. decimate again to pick one point per cubic 1m (for instance)
  // 3. now match the closest X points in 1 to those in 2, and generate surfel per point in 2.
  // 4. repeatedly: match surfels in cloud0 to those in cloud1
  // 5. solve a weighted least squares to find the transform of best fit
  // 6. apply the transform to cloud0 and save it out.

  vector<Matrix3d> matrices[2];

  vector<Vector3d> centroids[2];
  vector<Vector3d> normals[2];
  vector<Vector3d> widths[2];
  vector<bool> isPlane[2];
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
    isPlane[c].reserve(qSize);    
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

      Vector3d centroid;
      Vector3d width;
      Matrix3d mat;
      getSurfel(decimatedPoints, ids, centroid, width, mat);
      double q1 = width[0]/width[1];
      double q2 = width[1]/width[2];
      if (q2 < q1) // cylindrical
      {
        if (q2 > 0.5) // not cylinderical enough
          continue;
        // register two ellipsoids as the normal is ambiguous
        matrices[c].push_back(mat);
        centroids[c].push_back(centroid);
        normals[c].push_back(mat.col(2));       
        widths[c].push_back(width);  
        isPlane[c].push_back(false);
        if (c==1)
        {
          matrices[c].push_back(mat);
          centroids[c].push_back(centroid);
          normals[c].push_back(-mat.col(2));       
          widths[c].push_back(width);
          isPlane[c].push_back(false);  
        }
      }
      else
      {
        Vector3d normal = mat.col(0);
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
        getSurfel(decimatedPoints, ids, centroid, width, mat);
        normal = mat.col(0);
        double q1 = width[0]/width[1];
        
        if (q1 > 0.5) // not planar enough
          continue;
        if ((centroid - candidateStarts[i]).dot(normal) > 0.0)
          normal = -normal;     
        matrices[c].push_back(mat);
        centroids[c].push_back(centroid);
        normals[c].push_back(normal);       
        widths[c].push_back(width);  
        isPlane[c].push_back(true);
      }
    }
    draw.drawCloud(decimatedPoints, 0.5 + 0.4*(double)c, c);
  }
  draw.drawEllipsoids(centroids[1], matrices[1], widths[1], Vector3d(0,1,0), 1);

  // Now we have the ellipsoids, we need to do a nearest neighbour on this set of data, trying to match orientation as well as location
  struct Match
  {
    int ids[2];
    Vector3d normal;
  };
  vector<Match> matches;
  {
    // with these matches we now run the iterative reweighted least squares..
    Pose pose = Pose::identity();
    int maxIterations = 8;
    const double translationWeight = 0.4; // smaller finds matches further away
    const double maxNormalDifference = 0.5;
    for (int it = 0; it<maxIterations; it++)
    {
   //   if (it < 3) // the distribution of re-matching within the iterations is open to adjustment
      {
        vector<Vector3d> lineStarts;
        vector<Vector3d> lineEnds;
        int searchSize = 1;
        int qSize = centroids[0].size();
        int pSize = centroids[1].size();
        Nabo::NNSearchD *nns;
        MatrixXd pointsQ(7, qSize);
        for (int i = 0; i<qSize; i++)
        {
          Vector3d p = centroids[0][i]*translationWeight;
          p[2] *= 2.0; // doen't make much difference...
          pointsQ.col(i) << p, normals[0][i], isPlane[0][i]?1.0:0.0;
        }
        MatrixXd pointsP(7, pSize);
        for (int i = 0; i<pSize; i++)
        {
          Vector3d p = centroids[1][i]*translationWeight;
          p[2] *= 2.0;
          pointsP.col(i) << p, normals[1][i], isPlane[1][i]?1.0:0.0;
        }
        nns = Nabo::NNSearchD::createKDTreeLinearHeap(pointsP, 7);

        // Run the search
        MatrixXi indices;
        MatrixXd dists2;
        indices.resize(searchSize, qSize);
        dists2.resize(searchSize, qSize);
        nns->knn(pointsQ, indices, dists2, searchSize, 0.01, 0, maxNormalDifference);
        delete nns;
        matches.clear();

        for (int i = 0; i<qSize; i++)
        {
          // shall we pick only two-way matches? Not for now...
          for (int j = 0; j<searchSize && indices(j,i)>-1; j++)
          {
            Match match;
            match.ids[0] = i;
            match.ids[1] = indices(j,i);
            bool plane = isPlane[0][i];
            if (plane != isPlane[1][indices(j,i)])
              continue;
            Vector3d midNorm = (normals[0][i] + normals[1][indices(j,i)]).normalized();
            if (plane)
            {
              match.normal = midNorm;
              matches.push_back(match);
            }
            else
            {
              // a cylinder is like two normal constraints
              match.normal = midNorm.cross(Vector3d(1,2,3)).normalized();
              matches.push_back(match);
              match.normal = midNorm.cross(match.normal);
              matches.push_back(match);
            }
            lineStarts.push_back(centroids[0][i]);
            lineEnds.push_back(centroids[1][indices(j,i)]);
          }
        }
        draw.drawLines(lineStarts, lineEnds);
        draw.drawEllipsoids(centroids[0], matrices[0], widths[0], Vector3d(1,0,0), 0);
      }

      // don't go above 30*... or below 10*...
      double d = 20.0*(double)it/(double)maxIterations;
      const int stateSize = 12;
      Matrix<double, stateSize, stateSize> AtA;
      AtA.setZero();
      Matrix<double, stateSize, 1> AtB;
      AtB.setZero();
      double squareError = 0.0;
      for (int i = 0; i<(int)matches.size(); i++)
      {
        auto &match = matches[i];
        Vector3d pos[2] = {centroids[0][match.ids[0]], centroids[1][match.ids[1]]};
        double error = (pos[1] - pos[0]).dot(match.normal); // mahabolonis instead?
        double errorSqr;
        if (isPlane[0][match.ids[0]])
          errorSqr = sqr(error*translationWeight);
        else
        {
          Vector3d flat = pos[1] - pos[0];
          Vector3d norm = normals[0][match.ids[0]];
          flat -= norm * flat.dot(norm);
          errorSqr = (flat*translationWeight).squaredNorm();
        }
        // the normal difference is part of the error, 
        errorSqr += (normals[0][match.ids[0]] - normals[1][match.ids[1]]).squaredNorm();
        double weight = pow(max(1.0 - errorSqr/sqr(maxNormalDifference), 0.0), d*d);
        squareError += sqr(error);
        Matrix<double, stateSize, 1> At; // the Jacobian
        At.setZero();

        for (int i = 0; i<3; i++) // change in error with change in raycloud translation
          At[i] = match.normal[i];
        for (int i = 0; i<3; i++) // change in error with change in raycloud orientation
        {
          Vector3d axis(0,0,0);
          axis[i] = 1.0;
          At[3+i] = -(pos[0].cross(axis)).dot(match.normal);
        }
        if (!rigidOnly) // give the aligner a chance to rigidly align first
        {
          At[6] = sqr(pos[0][0]) * match.normal[0];
          At[7] = sqr(pos[0][0]) * match.normal[1];
          At[8] = sqr(pos[0][1]) * match.normal[0];
          At[9] = sqr(pos[0][1]) * match.normal[1];
          At[10] = pos[0][0]*pos[0][1] * match.normal[0];
          At[11] = pos[0][0]*pos[0][1] * match.normal[1];
        }
        AtA += At*weight*At.transpose();
        AtB += At*weight*error;
      }
      Matrix<double, stateSize, 1> x = AtA.ldlt().solve(AtB);
      cout << "rmse: " << sqrt(squareError/(double)matches.size()) << endl;
      cout << "least squares shift: " << x[0] << ", " << x[1] << ", " << x[2] << ", rotation: " << x[3] << ", " << x[4] << ", " << x[5] << endl;
      Vector3d rot(x[3], x[4], x[5]);
      Vector3d a(x[6],x[7],0), b(x[8],x[9],0), c(x[10],x[11],0);
      double angle = rot.norm();
      rot.normalize();
      Pose shift(Vector3d(x[0], x[1], x[2]), Quaterniond(AngleAxisd(angle, rot)));
      pose = shift * pose;
      for (int i = 0; i<(int)centroids[0].size(); i++)
      {
        Vector3d &pos = centroids[0][i];
        if (!rigidOnly)
          pos += a*sqr(pos[0]) + b*sqr(pos[1]) + c*pos[0]*pos[1];
        pos = shift*pos;
        normals[0][i] = shift.rotation * normals[0][i];
      }
      Quaterniond halfRot(AngleAxisd(angle/2.0, rot));
      for (auto &match: matches)
        match.normal = halfRot * match.normal;
      
      // TODO: transforming the whole cloud each time is a bit slow,
      // we should be able to concatenate these transforms and only apply them once at the end
      for (auto &end: aligner.clouds[0].ends)
      {
        if (!rigidOnly)
          end += a*sqr(end[0]) + b*sqr(end[1]) + c*end[0]*end[1];
        end = shift*end;      
      }
    }
  }

  aligner.clouds[0].save(fileStub + "_aligned.ply");  
  return true;
}
