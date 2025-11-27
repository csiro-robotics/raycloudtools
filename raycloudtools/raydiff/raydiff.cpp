// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <chrono>
#include <nabo/nabo.h>

#include "raylib/raycloud.h"
#include "raylib/rayparse.h"
#include "raylib/raycuboid.h"
#include "raylib/rayply.h"
// we look for a uniform distribution of best fit given a correlation dimension k
//  if the points are in a line   then k=1
//  if the points are in a plane  then k=2
//  if the points are in a volume then k=3
// this macro finds the k that best fit the uniform distribution. Otherwise we default to k=2
#define FIND_CORRELATION_DIMENSION 

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Difference between two ray clouds, differences coloured red, and similarity printed to screen. Optional visualisation." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydiff cloud1.ply cloud2.ply" << std::endl;
  std::cout << "                              --distance 0 - optional threshold in m for colouring differences. Default auto-detects distribution shoulder" << std::endl;
  std::cout << "                              --visualise  - open in the default visualisation tool" << std::endl;
  // clang-format on
  exit(exit_code);
}

void calcNearestNeighbourDistances(const ray::Cloud &cloud1, const ray::Cloud &cloud2, std::vector<double> &distances)
{
  Nabo::NNSearchD *nns;
 // Nabo::Parameters params("bucketSize", 8);
  int num_bounded = 0, num_bounded2 = 0;
  for (unsigned int i = 0; i < cloud1.ends.size(); i++)
    if (cloud1.rayBounded(i))
      num_bounded++;
  for (unsigned int i = 0; i < cloud2.ends.size(); i++)
    if (cloud2.rayBounded(i))
      num_bounded2++;
  Eigen::MatrixXd points_p(3, num_bounded);
  int j = 0;
  for (unsigned int i = 0; i < cloud1.ends.size(); i++) 
  {
    if (cloud1.rayBounded(i))
      points_p.col(j++) = cloud1.ends[i];
  }
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  Eigen::MatrixXd points_q(3, num_bounded2);
  j = 0;
  for (unsigned int i = 0; i < cloud2.ends.size(); i++) 
  {
    if (cloud2.rayBounded(i))
      points_q.col(j++) = cloud2.ends[i];
  }
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  // Run the search
  const int search_size = 1;
  indices.resize(search_size, num_bounded2);
  dists2.resize(search_size, num_bounded2);
  nns->knn(points_q, indices, dists2, search_size, ray::kNearestNeighbourEpsilon, Nabo::NNSearchD::SearchOptionFlags::SORT_RESULTS | Nabo::NNSearchD::SearchOptionFlags::ALLOW_SELF_MATCH);
  distances.resize(num_bounded2);
  for (int i = 0; i<num_bounded2; i++)
  {
    distances[i] = std::sqrt(dists2(0,i));
  }
  delete nns;
}

double getShoulder(double k, std::vector<double> sorted_dists, double &min_error_dist, double &similarity)
{
  int num = (int)sorted_dists.size();

  // 1. transform the result to make the distances uniform if their 3D points are uniformly distributed
  for (int i = 0; i<num; i++)
  {
    sorted_dists[i] = std::pow(sorted_dists[i], k); 
  }

  // 2. accumulate outside term backwards
  std::vector<double> outside_const(num+1, 0.0);
  std::vector<double> outside_linear(num+1, 0.0);
  std::vector<double> outside_square(num+1, 0.0);
  for (int i = num-1; i>=0; i--)
  {
    double d = sorted_dists.back() - sorted_dists[i];
    double yN = (num-1);
    double I = (double)i - yN;
    outside_const[i] = outside_const[i+1] + I*I;
    outside_linear[i] = outside_linear[i+1] + 2.0*I*d;
    outside_square[i] = outside_square[i+1] + d*d;      
  }

  // 3. accumulate linear term forwards, but store only best results
  double inside_const = 0.0;
  double inside_linear = 0.0;
  double inside_square = 0.0;
  double min_error_sqr = 0.0;
  double min_error_i = 0.0; 
  min_error_dist = 0.0;
  for (int i = 0; i<num; i++)
  {
    inside_const += ray::sqr((double)i);
    inside_linear -= 2.0*(double)i*sorted_dists[i];
    inside_square += ray::sqr(sorted_dists[i]);

    // ai^2 + bi + c = 0
    double yN = (double)(num-1);
    double d = sorted_dists.back() - sorted_dists[i];
    double square = outside_square[i]/ray::sqr(d);
    double linear = outside_linear[i]/d;
    double a = square;
    double b = -linear - 2.0*yN*square;
    double c = outside_const[i] + linear*yN + square*yN*yN; 

    // add the inside part
    a += inside_square/ray::sqr(sorted_dists[i]); // the division is to match the gradient to y
    b += inside_linear/sorted_dists[i];
    c += inside_const;

    double min_i = -b/(2.0*a);  // the height y of the integrated readings
    double error_sqr = a*min_i*min_i + b*min_i + c;
    if (error_sqr < min_error_sqr || i==0)
    {
      min_error_sqr = error_sqr;           // total square error in cumulative index, so measured in index offsets squared
      min_error_i = i;                     // number of points within uniform distribution of best fit
      min_error_dist = sorted_dists[i]; // the distance associated with the last index within the uniform distribution
    }
  }
  // 4. untransform the result:
  min_error_dist = std::pow(min_error_dist, 1.0/k);
  similarity = 100.0 * (double)min_error_i / (double)num;
  return std::sqrt(min_error_sqr);
}

int rayDiff(int argc, char *argv[])
{
  ray::FileArgument cloud1_name, cloud2_name;
  ray::OptionalFlagArgument visualise("visualise", 'v');
  ray::DoubleArgument distance_threshold(0.0, 1000.0, 0.0);
  ray::OptionalKeyValueArgument distance_option("distance", 'd', &distance_threshold);
  if (!ray::parseCommandLine(argc, argv, { &cloud1_name, &cloud2_name }, { &distance_option, &visualise }))
  {
    usage();
  }

  ray::Cloud cloud1, cloud2;
  if (!cloud1.load(cloud1_name.name()) || !cloud2.load(cloud2_name.name()))
  {
    usage();
  }

  // Primary goal - get an overall similarity percentage and return it as the return value, so it can be used in bash scripts
  // this needs to include colour and location, but be insensitive to changes in density. i.e. it needs to work usefully on a repeat scan.
  std::vector<double> dists_to_cloud2, dists_to_cloud1;
  std::cout << "calculating cloud2 neighbours of cloud1..." << std::endl;
  calcNearestNeighbourDistances(cloud1, cloud2, dists_to_cloud2);
  std::cout << "calculating cloud1 neighbours of cloud2..." << std::endl;
  calcNearestNeighbourDistances(cloud2, cloud1, dists_to_cloud1);

  // 1. max distance
  double max_dist = 0.0;
  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud2[i]);
  for (int i = 0; i<(int)dists_to_cloud1.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud1[i]);

  std::vector<double> sorted_dists = dists_to_cloud2;
  sorted_dists.insert(sorted_dists.end(), dists_to_cloud1.begin(), dists_to_cloud1.end());

  double min_error_dist = distance_threshold.value();
  double similarity = 0;

  if (distance_option.isSet())
  {
    std::cout << std::endl;
    std::cout << "Differences:" << std::setprecision(3) << std::fixed << std::endl;
    std::cout << std::endl;
    int inside_count = 0;
    for (int i = 0; i<(int)sorted_dists.size(); i++)
    {
      if (sorted_dists[i] <= min_error_dist)
      {
        inside_count++;
      }
    }
    similarity = 100.0*(double)inside_count / (int)sorted_dists.size();
    std::cout << " specified difference:    " << min_error_dist << " m,\t" << similarity << "% inside" << std::endl;
  }
  else
  {
    std::cout << "sorting differences..." << std::endl;
    std::sort(sorted_dists.begin(), sorted_dists.end());
    double min_error = 0;
    int min_i = 0;
    double k = 2.0; // good default as point clouds are typically surfaces
  #if defined FIND_CORRELATION_DIMENSION
    double errors[5];
    std::cout << "(uniform deviations: ";
    for (int i = 0; i<5; i++)
    {
      double k = 1.0 + (double)i/2.0;
      errors[i] = getShoulder(k, sorted_dists, min_error_dist, similarity);
      std::cout << "k " << k << ": " << errors[i];
      if (errors[i] < min_error || i==0)
      {
        min_error = errors[i];
        min_i = i;
      }
    }
    std::cout << ")" << std::endl;
    k = 1.0 + min_i/2.0;
    min_i = std::max(1, std::min(min_i, 3)); // so we can interpolate

    double y0 = errors[min_i-1];
    double y1 = errors[min_i];
    double y2 = errors[min_i+1];
    double den = y2 - 2.0*y1 + y0;
    if (den > 0.0)
    {
      double xmin = (y2 - 4.0*y1 + 3.0*y0)/(2.0*den);
      xmin = std::max(0.0, std::min(xmin, 2.0)); // clamp
      double I = (double)min_i - 1.0 + xmin;
      k = 1.0 + I/2.0;
    }
  #endif

    std::cout << std::endl;
    std::cout << "Differences:" << std::setprecision(3) << std::fixed << std::endl;
    std::cout << std::endl;
    min_error = getShoulder(k, sorted_dists, min_error_dist, similarity);
    std::cout << " shoulder difference:     " << min_error_dist << " m,\t" << similarity << "% inside, correlation dimension: " << k << std::endl;
  }
  std::cout << " median difference:       " << sorted_dists[sorted_dists.size()/2] << " m" << std::endl;
  std::cout << " max difference:          " << max_dist << " m" << std::endl;
  std::cout << std::endl;

  std::cout << "saving out coloured red for matches beyond shoulder difference:" << std::endl;
  // now render visuals
  int j = 0;
  Eigen::Vector3d diff_col(255,0,0);
  for (int i = 0; i<(int)cloud1.ends.size(); i++)
  {
    if (cloud1.rayBounded(i))
    {
      if (dists_to_cloud1[j] > min_error_dist)
      {
        double proximity = min_error_dist / dists_to_cloud1[j];
        Eigen::Vector3d col(cloud1.colours[i].red, cloud1.colours[i].green, cloud1.colours[i].blue);
        Eigen::Vector3d new_col = diff_col + (col - diff_col) * proximity;
        cloud1.colours[i] = ray::RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), cloud1.colours[i].alpha);
      }
      j++;
    }
  }
  cloud1.save(cloud1_name.nameStub() + "_diff.ply");
  j = 0;
  for (int i = 0; i<(int)cloud2.ends.size(); i++)
  {
    if (cloud2.rayBounded(i))
    {
      if (dists_to_cloud2[j] > min_error_dist)
      {
        double proximity = min_error_dist / dists_to_cloud2[j];
        Eigen::Vector3d col(cloud2.colours[i].red, cloud2.colours[i].green, cloud2.colours[i].blue);
        Eigen::Vector3d new_col = diff_col + (col - diff_col) * proximity;
        cloud2.colours[i] = ray::RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), cloud2.colours[i].alpha);
      }
      j++;
    }
  }
  cloud2.save(cloud2_name.nameStub() + "_diff.ply");

  if (visualise.isSet())
  {
    std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_name.nameStub() + "_diff.ply " + cloud2_name.nameStub() + "_diff.ply";
    system(command.c_str());  
  }
  return (int)similarity;
  return 1;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDiff, argc, argv);
}