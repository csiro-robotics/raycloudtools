// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayaxisalign.h"
#include "rayutils.h"
#include "raycloud.h"
#include "raycloudwriter.h"
#include "rayunused.h"
#include "rayply.h"

namespace ray
{
namespace
{
/// These constants are the angular and amplitude (distance from centre) resolutions for the radon transform. 
/// Too small and the alignment is imprecise, too large and it fails to align on noisy planes (such as a crop row).
/// We may ultimately make these a function of the cloud's size or density, but for this iteration of the axis alignment function, they are fixed.
const int ang_res = 256; // ang_res must be divisible by 2
const int amp_res = 256;  
}

/// the output file @c out_file is the input file @c in_file, transformed by @c pose
bool transformCloud(const std::string &in_file, const std::string &out_file, const Pose &pose)
{
  CloudWriter writer;
  if (!writer.begin(out_file))
    return false;
  Cloud chunk;

  auto transform = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours)
  {
    chunk.clear();
    chunk.times = times;
    chunk.colours = colours;
    for (size_t i = 0; i<ends.size(); i++)
    {
      chunk.starts.push_back(pose * starts[i]);
      chunk.ends.push_back(pose * ends[i]);
    }
    writer.writeChunk(chunk);
  };

  if (!Cloud::read(in_file, transform))
    return false;

  writer.end();
  return true;
}

// just a quadratic maximum -b/2a for heights y0,y1,y2
double peak(double y0, double y1, double y2)
{
  return 0.5 + 0.5 * (y0 - y2) / (y0 + y2 - 2.0 * y1);  
}

Pose estimatePose(const Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> &position_accumulator, const Cloud::Info &info)
{
  const Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  const Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;
  const Eigen::Vector3d mid_bound = (min_bound + max_bound)/2.0;

  // radius used for maximum amplitude (distance from mid_bound)
  const double radius = 0.5 * std::sqrt(ray::sqr(max_bound[0]-min_bound[0]) + ray::sqr(max_bound[1]-min_bound[1]));
  const double eps = 0.0001; // avoids the most distant point 
  Eigen::ArrayXXd weights(ang_res, amp_res);
  weights.fill(0);
  // process the data into a Euclidean transform called pose:
  for (int ii = 0; ii<amp_res; ii++)
  {
    for (int jj = 0; jj<amp_res; jj++)
    {
      const Eigen::Vector3d &accumulator = position_accumulator(ii, jj);
      double w = accumulator[2];
      if (w == 0.0)
        continue;
      Eigen::Vector2d centroid(accumulator[0]/w, accumulator[1]/w);
      double angle = atan2(centroid[0], centroid[1]);
      double amplitude = std::sqrt(centroid[0]*centroid[0] + centroid[1]*centroid[1]) / radius;

      // now draw the sine wave for this point. (should it be anti-aliased?)
      for (int i = 0; i<ang_res; i++)
      {
        double ang = kPi * (double)i/(double)ang_res;
        double y = amplitude * std::sin(ang + angle);
        double x = ((double)amp_res-1.0-eps) * (0.5 + 0.5*y);
        int j = (int)x;
        double blend = x - (double)j;
        
        weights(i, j) += (1.0-blend)*w;
        weights(i, j+1) += blend*w;
      }
    } 
  }
  // now find heighest weight:    
  int max_i = 0, max_j = 0;
  double max_weight = 0.0;
  for (int i = 0; i<ang_res; i++)
  {
    for (int j = 0; j<amp_res; j++)
    {
      if (weights(i, j) > max_weight)
      {
        max_i = i;
        max_j = j;
        max_weight = weights(i, j);
      }
    }
  }
  std::cout << "max_i: " << max_i << ", max_j: " << max_j << std::endl;

  // now interpolate the location
  double angle = (double)max_i + peak(weights((max_i + ang_res-1)%ang_res, max_j), weights(max_i, max_j), weights((max_i+1)%ang_res, max_j));
  angle *= kPi / (double)ang_res;
  std::cout << "principle angle: " << angle << ", direction vec: " << std::cos(angle) << ", " << std::sin(angle) << std::endl;

  double amp = (double)max_j + peak(weights(max_i, std::max(0, max_j-1)), weights(max_i, max_j), weights(max_i, std::min(max_j+1, amp_res-1)));
  amp = radius * ((2.0 * amp/(double)amp_res) - 1.0);
  std::cout << "distance along direction: " << amp << ", radius: " << radius << "maxj/ampres: " << max_j/(double)amp_res << std::endl;

  Eigen::Vector2d line_vector = amp * Eigen::Vector2d(std::cos(angle), std::sin(angle));

  // now find the orthogonal best edge. Simply the greatest weight along the orthogonal angle
  int orth_i = (max_i + ang_res/2)%ang_res;
  double max_w = 0.0;
  int max_orth_j = 0;
  for (int j = 0; j<amp_res; j++)
  {
    if (weights(orth_i, j) > max_w)
    {
      max_w = weights(orth_i, j);
      max_orth_j = j;
    }
  }
  // now interpolate this orthogonal direction:
  double amp2 = (double)max_orth_j + peak(weights(orth_i, std::max(0, max_orth_j-1)), weights(orth_i, max_orth_j), weights(orth_i, std::min(max_orth_j+1, amp_res-1)));
  amp2 = radius * ((2.0 * amp2/(double)amp_res) - 1.0);
  std::cout << "orthi: " << orth_i << ", max_i: " << max_i << std::endl;
  if (orth_i < max_i) // the %res earlier puts it in antiphase (since we only have 180 degrees per weights map)
    amp2 = -amp2;     // so negate the amplitude here

  Eigen::Vector2d line2_vector = amp2 * Eigen::Vector2d(std::cos(angle+kPi/2.0), std::sin(angle+kPi/2.0));
  std::cout << "line2_vector: " << line2_vector.transpose() << std::endl;
  std::cout << "amp: " << amp << ", amp2: " << amp2 << std::endl;

  Eigen::Vector2d centre = line_vector + line2_vector;
  Eigen::Quaterniond id(1,0,0,0);
  Pose to_mid(-mid_bound - Eigen::Vector3d(centre[0], centre[1], 0), id);
  Pose rot(Eigen::Vector3d::Zero(), Eigen::Quaterniond(Eigen::AngleAxisd(-angle, Eigen::Vector3d(0,0,1))));
  Pose pose = rot * to_mid;

  // now rotate the cloud into the positive x/y axis, depending on which is the stronger signal
  Eigen::Vector3d mid_point = pose * info.centroid;
  if (std::abs(mid_point[0]) > std::abs(mid_point[1]))
  {
    if (mid_point[0] < 0.0)
      pose = Pose(Eigen::Vector3d(0,0,0), Eigen::Quaterniond(0,0,0,1)) * pose; // 180 degree yaw
  }
  else
  {
    if (mid_point[1] < 0.0)
      pose = Pose(Eigen::Vector3d(0,0,0), Eigen::Quaterniond(0,0,0,1)) * pose; // 180 degree yaw
  }
  return pose;
}

// A radon transform is used as follows:
// 1. we quantise the cloud into a 2D grid of centroids, weighted by number of end points within the cell
// 2. for each cell, render a sine wave to weights image
// 3. find the maximum weight in the image (the new y axis) and interpolate the angle and position using its neighbours
// 4. find the highest orthogonal weight (the new x axis) and interpolate its position
// 5. quantise density vertically into an array to get the strongest ground height signal, interpolating the max value
// The advantage of the radon transform is that it does not require normals (unreliable on vegetation), it can work on 
// fairly noisy planes (such as a vineyard row), and it should parallelise well.
bool alignCloudToAxes(const std::string &cloud_name, const std::string &aligned_file)
{
  // Calculate extents:
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
    return false;
  const Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  const Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;
  const Eigen::Vector3d mid_bound = (min_bound + max_bound)/2.0;

  // fill in the arrays
  double eps = 0.0001;

  // first convert the cloud into a weighted centroid field. I'm using element [2] for the weight.
  Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> position_accumulator(amp_res, amp_res);
  position_accumulator.fill(Eigen::Vector3d::Zero());
  double step_x = ((double)amp_res-1.0-eps) / (max_bound[0]-min_bound[0]);
  double step_y = ((double)amp_res-1.0-eps) / (max_bound[1]-min_bound[1]);
  // and set up the vertical arrays too
  double vertical_weights[amp_res]; // TODO: height is usually much less than width... more constant voxel size?
  std::memset(vertical_weights, 0, sizeof(double)*amp_res);
  double step_z = ((double)amp_res - 1.0 - eps) / (max_bound[2] - min_bound[2]);

  auto fill_arrays = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    RAYLIB_UNUSED(starts);
    RAYLIB_UNUSED(times);
    RAYLIB_UNUSED(colours);
    for (size_t e = 0; e<ends.size(); e++)
    {
      if (colours[e].alpha == 0) // unbounded
        continue;
      Eigen::Vector3d index = ends[e] - min_bound;
      index[0] *= step_x;
      index[1] *= step_y;
      index[2] *= step_z;
      Eigen::Vector3d pos = ends[e] - mid_bound;
      pos[2] = 1.0;
      position_accumulator((int)index[0], (int)index[1]) += pos;

      int k = (int)index[2];
      double blend = index[2] - (double)k;
      vertical_weights[k] += (1.0 - blend);
      vertical_weights[k+1] += blend;      
    }
  };
  if (!Cloud::read(cloud_name, fill_arrays))
    return false;

  Pose pose = estimatePose(position_accumulator, info);

  // now get the vertical displacement
  int max_k = 0;
  double max_wt = 0;
  for (int k = 0; k<amp_res; k++)
  {
    if (vertical_weights[k] > max_wt) 
    {
      max_wt = vertical_weights[k];
      max_k = k;
    }
  }
  double k2 = (double)max_k + peak(vertical_weights[std::max(0, max_k-1)], vertical_weights[max_k], vertical_weights[std::min(max_k+1, amp_res-1)]);
  double height = (k2 / step_z) + min_bound[2];
  pose.position[2] = -height;

  std::cout << "pose: " << pose.position.transpose() << ", q: " << pose.rotation.w() << ", " << pose.rotation.x() << ", " << pose.rotation.y() << ", " << pose.rotation.z() << std::endl;

  // finally, transform the cloud:
  return transformCloud(cloud_name, aligned_file, pose);
}

} // namespace ray
