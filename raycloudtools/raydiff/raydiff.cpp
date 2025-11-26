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

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Difference between two ray clouds, differences coloured red, and similarity printed to screen. Optional visualisation." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydiff cloud1.ply cloud2.ply" << std::endl;
  std::cout << "                              --visualise     - open in the default visualisation tool" << std::endl;
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

int rayDiff(int argc, char *argv[])
{
  ray::FileArgument cloud1_name, cloud2_name;
  ray::OptionalFlagArgument visualise("--visualise", 'v');
  ray::IntArgument type(0, 10, 0);
  ray::OptionalKeyValueArgument type_option("--type", 't', &type);
  if (!ray::parseCommandLine(argc, argv, { &cloud1_name, &cloud2_name }, { &type_option, &visualise }))
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
  calcNearestNeighbourDistances(cloud1, cloud2, dists_to_cloud2);
  calcNearestNeighbourDistances(cloud2, cloud1, dists_to_cloud1);

  // 1. max distance
  double max_dist = 0.0;
  // 2. percentage of cloud within tolerance

  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud2[i]);
  for (int i = 0; i<(int)dists_to_cloud1.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud1[i]);
  int num1 = (int)dists_to_cloud1.size();
  int num2 = (int)dists_to_cloud2.size();

  // 3. median distance
  std::vector<double> sorted_dists = dists_to_cloud2;
  sorted_dists.insert(sorted_dists.end(), dists_to_cloud1.begin(), dists_to_cloud1.end());
  double median_dist = ray::median(sorted_dists);

  std::sort(sorted_dists.begin(), sorted_dists.end());
  int num = (int)sorted_dists.size();
  // 1. transform the result to make the distances uniform if their 3D points are uniformly distributed
  for (int i = 0; i<num; i++)
  {
    double d = sorted_dists[i];
    sorted_dists[i] = d*d*d; // would use d*d if the points were planar
  }
 
  bool piecewise_linear = type.value() == 1;

  // 2. accumulate uniform term backwards
  std::vector<double> uniform_const(num+1, 0.0);
  std::vector<double> uniform_linear(num+1, 0.0);
  std::vector<double> uniform_square(num+1, 0.0);
  for (int i = num-1; i>=0; i--)
  {
    if (piecewise_linear)
    {
      double d = sorted_dists.back() - sorted_dists[i];
      double yN = (num-1);
      double I = (double)i - yN;
      uniform_const[i] = uniform_const[i+1] + I*I;
      uniform_linear[i] = uniform_linear[i+1] + 2.0*I*d;
      uniform_square[i] = uniform_square[i+1] + d*d;      
    }
    else
    {
      uniform_const[i] = uniform_const[i+1] + ray::sqr((double)i);
      uniform_linear[i] = uniform_linear[i+1] - 2.0*(double)i;
      uniform_square[i] = uniform_square[i+1] + 1.0;
    }
  }

  // 3. accumulate linear term forwards, but store only best results
  double linear_const = 0.0;
  double linear_linear = 0.0;
  double linear_square = 0.0;
  double min_error_sqr = 0.0;
  double min_error_i = 0.0; 
  double min_error_dist = 0.0;
  for (int i = 0; i<num; i++)
  {
    linear_const += ray::sqr((double)i);
    linear_linear -= 2.0*(double)i*sorted_dists[i];
    linear_square += ray::sqr(sorted_dists[i]);

    // ai^2 + bi + c = 0
    double yN = (double)(num-1);
    double a = 0, b = 0, c = 0;
    if (piecewise_linear)
    {
      double d = sorted_dists.back() - sorted_dists[i];

      double square = uniform_square[i]/ray::sqr(d);
      double linear = uniform_linear[i]/d;
      a = square;
      b = -linear - 2.0*yN*square;
      c = uniform_const[i] + linear*yN + square*yN*yN; 

      // add the inside part
      a += linear_square/ray::sqr(sorted_dists[i]); // the division is to match the gradient to y
      b += linear_linear/sorted_dists[i];
      c += linear_const;
    }
    else
    {
      a = uniform_square[i] + linear_square/ray::sqr(sorted_dists[i]); // the division is to match the gradient to y
      b = uniform_linear[i] + linear_linear/sorted_dists[i];
      c = uniform_const[i] + linear_const;
    }

    double min_i = -b/(2.0*a);
    double error_sqr = a*min_i*min_i + b*min_i + c;
/*
    double grad = min_i / sorted_dists[i];
    double error_sqr1 = linear_square*grad*grad + linear_linear*grad + linear_const;
    double grad2 = (yN - min_i) / (sorted_dists.back() - sorted_dists[i]);
    double error_sqr2 = uniform_square[i]*grad2*grad2 + uniform_linear[i]*grad2 + uniform_const[i];
    double sum = error_sqr1 + error_sqr2;
    if (std::abs(sum/error_sqr - 1.0) > 0.01)
    {
      std::cout << "these shouldn't be different" << std::endl;
    }*/

    if (error_sqr < min_error_sqr || i==0)
    {
      min_error_sqr = error_sqr;           // total square error in cumulative index, so measured in index offsets squared
      min_error_i = i;                     // number of points within uniform distribution of best fit
      min_error_dist = sorted_dists[i]; // the distance associated with the last index within the uniform distribution
    }
  }
  // 4. untransform the result:
  min_error_dist = std::pow(min_error_dist, 1.0/3.0);

  // print results...
  double similarity = 100.0 * (double)min_error_i / (double)num;
  std::cout << "shoulder distance: " << min_error_dist << " m, within: " << similarity << "%" << std::endl;
  std::cout << "median difference: " << median_dist << " m, max difference: " << max_dist << " m" << std::endl;

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
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDiff, argc, argv);
}