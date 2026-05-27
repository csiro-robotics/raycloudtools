// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tom Lowe
#include <nabo/nabo.h>
#include "raydiffer.h"
#include "raylib/raycloud.h"
#include "raylib/rayply.h"
#include "raylib/raycloudwriter.h"
// we look for a uniform distribution of best fit given a correlation dimension k
//  if the points are in a line   then k=1
//  if the points are in a plane  then k=2
//  if the points are in a volume then k=3
// this macro finds the k that best fit the uniform distribution. Otherwise we default to k=2
#define FIND_CORRELATION_DIMENSION 

namespace ray
{
void calcNearestNeighbourDistances(const std::vector<Eigen::Vector3f> &cloud1, const std::vector<Eigen::Vector3f> &cloud2, std::vector<float> &distances)
{
  Nabo::NNSearchF *nns;
 // Nabo::Parameters params("bucketSize", 8);
  int num_bounded = 0, num_bounded2 = 0;
  for (unsigned int i = 0; i < cloud1.size(); i++)
    if (cloud1[i][0] != std::numeric_limits<float>::max())
      num_bounded++;
  for (unsigned int i = 0; i < cloud2.size(); i++)
    if (cloud2[i][0] != std::numeric_limits<float>::max())
      num_bounded2++;
  Eigen::MatrixXf points_p(3, num_bounded);
  int j = 0;
  for (unsigned int i = 0; i < cloud1.size(); i++) 
  {
    if (cloud1[i][0] != std::numeric_limits<float>::max())
      points_p.col(j++) = cloud1[i];
  }
  nns = Nabo::NNSearchF::createKDTreeLinearHeap(points_p, 3);

  Eigen::MatrixXf points_q(3, num_bounded2);
  j = 0;
  for (unsigned int i = 0; i < cloud2.size(); i++) 
  {
    if (cloud2[i][0] != std::numeric_limits<float>::max())
      points_q.col(j++) = cloud2[i];
  }
  Eigen::MatrixXi indices;
  Eigen::MatrixXf dists2;
  // Run the search
  const int search_size = 1;
  indices.resize(search_size, num_bounded2);
  dists2.resize(search_size, num_bounded2);
  nns->knn(points_q, indices, dists2, search_size, (float)kNearestNeighbourEpsilon, Nabo::NNSearchF::SearchOptionFlags::SORT_RESULTS | Nabo::NNSearchF::SearchOptionFlags::ALLOW_SELF_MATCH);
  delete nns;
  distances.resize(num_bounded2);
  for (int i = 0; i<num_bounded2; i++)
  {
    distances[i] = (float)std::sqrt(dists2(0,i));
  }
}

void getDistancesBetweenPoints(const std::vector<Eigen::Vector3f> &points1, const std::vector<Eigen::Vector3f> &points2, std::vector<float> &dists_to_cloud1, std::vector<float> &dists_to_cloud2)
{
  std::cout << "calculating cloud2 neighbours of cloud1..." << std::endl;
  calcNearestNeighbourDistances(points1, points2, dists_to_cloud2);
  std::cout << "calculating cloud1 neighbours of cloud2..." << std::endl;
  calcNearestNeighbourDistances(points2, points1, dists_to_cloud1);
}

double getShoulder(double k, std::vector<float> sorted_dists, double &min_error_dist, double &similarity)
{
  int num = (int)sorted_dists.size();

  // 1. transform the result to make the distances uniform if their 3D points are uniformly distributed
  for (int i = 0; i<num; i++)
  {
    sorted_dists[i] = (float)std::pow((double)sorted_dists[i], k); 
  }

  // 2. accumulate outside term backwards
  std::vector<float> outside_const(num+1, 0.0);
  std::vector<float> outside_linear(num+1, 0.0);
  std::vector<float> outside_square(num+1, 0.0);
  for (int i = num-1; i>=0; i--)
  {
    double d = sorted_dists.back() - sorted_dists[i];
    double yN = (num-1);
    double I = (float)i - yN;
    outside_const[i] = (float)((double)outside_const[i+1] + I*I);
    outside_linear[i] = (float)((double)outside_linear[i+1] + 2.0*I*d);
    outside_square[i] =(float)((double)outside_square[i+1] + d*d);      
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
    inside_const += sqr((double)i);
    inside_linear -= 2.0*(double)i*sorted_dists[i];
    inside_square += sqr(sorted_dists[i]);

    // ai^2 + bi + c = 0
    double yN = (double)(num-1);
    double d = sorted_dists.back() - sorted_dists[i];
    double square = outside_square[i]/sqr(d);
    double linear = outside_linear[i]/d;
    double a = square;
    double b = -linear - 2.0*yN*square;
    double c = outside_const[i] + linear*yN + square*yN*yN; 

    // add the inside part
    a += inside_square/sqr(sorted_dists[i]); // the division is to match the gradient to y
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

double printDistanceStatistics(const std::vector<float> &dists_to_cloud1, const std::vector<float> &dists_to_cloud2, double &distance_threshold)
{
  // 1. max distance
  float max_dist = 0.0;
  for (int i = 0; i<(int)dists_to_cloud2.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud2[i]);
  for (int i = 0; i<(int)dists_to_cloud1.size(); i++)
    max_dist = std::max(max_dist, dists_to_cloud1[i]);

  std::vector<float> sorted_dists = dists_to_cloud2;
  sorted_dists.insert(sorted_dists.end(), dists_to_cloud1.begin(), dists_to_cloud1.end());

  double similarity = 0;

  if (distance_threshold > 0.0)
  {
    std::cout << std::endl;
    std::cout << "Differences:" << std::setprecision(3) << std::fixed << std::endl;
    std::cout << std::endl;
    int inside_count = 0;
    for (int i = 0; i<(int)sorted_dists.size(); i++)
    {
      if (sorted_dists[i] <= distance_threshold)
      {
        inside_count++;
      }
    }
    similarity = 100.0*(double)inside_count / (int)sorted_dists.size();
    std::cout << " specified difference:    " << distance_threshold << " m,\t" << similarity << "% inside" << std::endl;
  }
  else
  {
    // prune out distance 0
    for (int i = (int)sorted_dists.size()-1; i>=0; i--)
    {
      if (sorted_dists[i] <= 0.0)
      {
        sorted_dists[i] = sorted_dists.back();
        sorted_dists.pop_back();
      }
    }
    if (sorted_dists.empty())
    {
      std::cerr << "No difference detect between files. Exiting" << std::endl;
      return 100.0;
    }
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
      errors[i] = getShoulder(k, sorted_dists, distance_threshold, similarity);
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
    min_error = getShoulder(k, sorted_dists, distance_threshold, similarity);
    std::cout << " shoulder difference:     " << distance_threshold << " m,\t" << similarity << "% inside, correlation dimension: " << k << std::endl;
  }
  std::cout << " median difference:       " << sorted_dists[sorted_dists.size()/2] << " m" << std::endl;
  std::cout << " max difference:          " << max_dist << " m" << std::endl;
  std::cout << std::endl;
  return similarity;
}

bool writeDifferencesToRayClouds(const std::string &cloud1_namestub, const std::string &cloud2_namestub, std::vector<float> &dists_to_cloud1, 
   std::vector<float> &dists_to_cloud2, double dist_threshold, bool individual_files, bool visualise)
{
  std::cout << "saving out differences, coloured red for differences to " << cloud1_namestub << ".ply and green for differences in " << cloud2_namestub << ".ply" << std::endl;

  // now render visuals
  CloudWriter writer;
  if (!writer.begin(cloud1_namestub + "_diff.ply"))
    return false;

  // By maintaining these buffers below, we avoid almost all memory fragmentation
  Cloud chunk;
  std::map<Eigen::Vector3i, Eigen::Vector2i, Vector3iLess> voxel_map;
  std::vector<Eigen::Vector3i> samples;

  int j = 0;
  Eigen::Vector3d diff_col(255,0,0);
  std::vector<float> *dists = &dists_to_cloud1;
  const float eps = 1e-8f; // in case there is inaccuracy in the KNN distance estimation for co-located point pairs
  auto colour = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                    std::vector<double> &times, std::vector<RGBA> &colours) 
  {
    // firstly we store a count per cell
    chunk.clear();
    for (size_t i = 0; i<ends.size(); i++)
    {
      if (colours[i].alpha > 0)
      {
        if ((*dists)[j] > eps)
        {
          chunk.addRay(starts[i], ends[i], times[i], colours[i]);
          if ((*dists)[j] > dist_threshold)
          {
            double proximity = dist_threshold / (*dists)[j];
            Eigen::Vector3d col(colours[i].red, colours[i].green, colours[i].blue);
            Eigen::Vector3d new_col = diff_col + (col - diff_col) * proximity;
            chunk.colours.back() = RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), colours[i].alpha);
          }
        }
        j++;
      }
      else
        chunk.addRay(starts[i], ends[i], times[i], colours[i]);
    }
    writer.writeChunk(chunk);
  };

  if (!Cloud::read(cloud1_namestub + ".ply", colour))
    return false;
  j = 0;
  dists_to_cloud1.clear();
  dists_to_cloud1.shrink_to_fit();
  Eigen::Vector3d diff2_col(0,255,0);

  dists = &dists_to_cloud2;
  if (!individual_files)
  {
    auto colour2 = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<RGBA> &colours) 
    {
      // firstly we store a count per cell
      chunk.clear();
      for (size_t i = 0; i<ends.size(); i++)
      {
        if (colours[i].alpha > 0)
        {
          if ((*dists)[j] > eps)
          {
            chunk.addRay(starts[i], ends[i], times[i], colours[i]);
            if ((*dists)[j] > dist_threshold)
            {
              double proximity = dist_threshold / (*dists)[j];
              Eigen::Vector3d col(colours[i].red, colours[i].green, colours[i].blue);
              Eigen::Vector3d new_col = diff2_col + (col - diff2_col) * proximity;
              chunk.colours[i] = RGBA((uint8_t)(new_col[0]+0.5), (uint8_t)(new_col[1]+0.5), (uint8_t)(new_col[2]+0.5), colours[i].alpha);
            }
          }
          j++;
        }
        else
          chunk.addRay(starts[i], ends[i], times[i], colours[i]);
      }
      writer.writeChunk(chunk);
    };

    if (!Cloud::read(cloud2_namestub + ".ply", colour2))
      return false;
    writer.end();

    if (visualise)
    {
      std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_namestub + "_diff.ply";
      system(command.c_str());  
    }
  }
  else
  {
    writer.end();
    diff_col = diff2_col;
    if (!writer.begin(cloud2_namestub + "_diff.ply"))
      return false;

    if (!Cloud::read(cloud2_namestub + ".ply", colour))
      return false;
    writer.end();
    if (visualise)
    {
      std::string command = std::string(VISUALISE_TOOL) + std::string(" ") + cloud1_namestub + "_diff.ply " + cloud2_namestub + "_diff.ply";
      system(command.c_str());  
    }
  }
  return true;
}

}  // namespace ray
