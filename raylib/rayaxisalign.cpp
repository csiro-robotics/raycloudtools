// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayaxisalign.h"
#include "rayutils.h"
#include "raycloud.h"
#include "rayunused.h"
#include "rayply.h"

namespace ray
{
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
  const Eigen::Vector3d mid = (min_bound + max_bound)/2.0;

  // fill in the arrays
  const int ang_res = 256, amp_res = 256; // ang_res must be divisible by 2
  double weights[ang_res][amp_res];
  std::memset(weights, 0, sizeof(double)*ang_res*amp_res);
  double radius = 0.5 * std::sqrt(2.0) * std::max(max_bound[0]-min_bound[0], max_bound[1]-min_bound[1]);
  double eps = 0.0001;

  // first convert the cloud into a weighted centroid field. I'm using element [2] for the weight.
  Eigen::Vector3d ps[amp_res][amp_res];
  std::memset(ps, 0, sizeof(Eigen::Vector3d)*amp_res*amp_res);
  double step_x = ((double)amp_res-1.0-eps) / (max_bound[0]-min_bound[0]);
  double step_y = ((double)amp_res-1.0-eps) / (max_bound[1]-min_bound[1]);
  // and set up the vertical arrays too
  double ws[amp_res]; // TODO: height is usually much less than width... more constant voxel size?
  std::memset(ws, 0, sizeof(double)*amp_res);
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
      Eigen::Vector3d pos = ends[e] - mid;
      pos[2] = 1.0;
      ps[(int)index[0]][(int)index[1]] += pos;

      int k = (int)index[2];
      double blend = index[2] - (double)k;
      ws[k] += (1.0 - blend);
      ws[k+1] += blend;      
    }
  };
  if (!Cloud::read(cloud_name, fill_arrays))
    return false;

  // process the data into a Euclidean transform called pose:
  for (int ii = 0; ii<amp_res; ii++)
  {
    for (int jj = 0; jj<amp_res; jj++)
    {
      Eigen::Vector3d &p = ps[ii][jj];
      double w = p[2];
      if (w == 0.0)
        continue;
      Eigen::Vector2d pos(p[0]/w, p[1]/w);
      double angle = atan2(pos[0], pos[1]);
      double amp = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]) / radius;

      // now draw the sine wave for this point. (should it be anti-aliased?)
      for (int i = 0; i<ang_res; i++)
      {
        double ang = kPi * (double)i/(double)ang_res;
        double y = amp * std::sin(ang + angle);
        double x = ((double)amp_res-1.0-eps) * (0.5 + 0.5*y);
        int j = (int)x;
        double blend = x - (double)j;
        
        weights[i][j] += (1.0-blend)*w;
        weights[i][j+1] += blend*w;
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
      if (weights[i][j] > max_weight)
      {
        max_i = i;
        max_j = j;
        max_weight = weights[i][j];
      }
    }
  }
  std::cout << "max_i: " << max_i << ", max_j: " << max_j << std::endl;

  // now interpolate the location
  double x0 = weights[(max_i + ang_res-1)%ang_res][max_j];
  double x1 = weights[max_i][max_j];
  double x2 = weights[(max_i+1)%ang_res][max_j];
  double angle = (double)max_i+0.5 + 0.5 * (x0 - x2) / (x0 + x2 - 2.0 * x1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
  angle *= kPi / (double)ang_res;
  std::cout << "principle angle: " << angle << ", direction vec: " << std::cos(angle) << ", " << std::sin(angle) << std::endl;

  double y0 = weights[max_i][std::max(0, max_j-1)];
  double y1 = weights[max_i][max_j];
  double y2 = weights[max_i][std::min(max_j+1, amp_res-1)];
  double amp = (double)max_j+0.5 + 0.5 * (y0 - y2) / (y0 + y2 - 2.0 * y1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
  amp = radius * ((2.0 * amp/(double)amp_res) - 1.0);
  std::cout << "distance along direction: " << amp << ", radius: " << radius << "maxj/ampres: " << max_j/(double)amp_res << std::endl;

  Eigen::Vector2d line_vector = amp * Eigen::Vector2d(std::cos(angle), std::sin(angle));

  // now find the orthogonal best edge. Simply the greatest weight along the orthogonal angle
  int orth_i = (max_i + ang_res/2)%ang_res;
  double max_w = 0.0;
  int max_orth_j = 0;
  for (int j = 0; j<amp_res; j++)
  {
    if (weights[orth_i][j] > max_w)
    {
      max_w = weights[orth_i][j];
      max_orth_j = j;
    }
  }
  // now interpolate this orthogonal direction:
  double z0 = weights[orth_i][std::max(0, max_orth_j-1)];
  double z1 = weights[orth_i][max_orth_j];
  double z2 = weights[orth_i][std::min(max_orth_j+1, amp_res-1)];
  double amp2 = (double)max_orth_j+0.5 + 0.5 * (z0 - z2) / (z0 + z2 - 2.0 * z1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
  amp2 = radius * ((2.0 * amp2/(double)amp_res) - 1.0);
  std::cout << "orthi: " << orth_i << ", max_i: " << max_i << std::endl;
  if (orth_i < max_i) // the %res earlier puts it in antiphase (since we only have 180 degrees per weights map)
    amp2 = -amp2;     // so negate the amplitude here

  Eigen::Vector2d line2_vector = amp2 * Eigen::Vector2d(std::cos(angle+kPi/2.0), std::sin(angle+kPi/2.0));
  std::cout << "line2_vector: " << line2_vector.transpose() << std::endl;
  std::cout << "amp: " << amp << ", amp2: " << amp2 << std::endl;

  // great, we now have the angle and cross-over position. Next is to flip it so the largest side is positive
  Eigen::Vector2d centre = line_vector + line2_vector;
  Eigen::Quaterniond id(1,0,0,0);
  Pose to_mid(-mid - Eigen::Vector3d(centre[0], centre[1], 0), id);
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

  // now get the vertical displacement
  int max_k = 0;
  double max_wt = 0;
  for (int k = 0; k<amp_res; k++)
  {
    if (ws[k] > max_wt) 
    {
      max_wt = ws[k];
      max_k = k;
    }
  }
  double w0 = ws[std::max(0, max_k-1)];
  double w1 = ws[max_k];
  double w2 = ws[std::min(max_k+1, amp_res-1)];
  double k2 = (double)max_k+0.5 + 0.5 * (w0 - w2) / (w0 + w2 - 2.0 * w1);  // just a quadratic maximum -b/2a for heights y0,y1,y2
  double height = (k2 / step_z) + min_bound[2];
  pose.position[2] = -height;

  std::cout << "pose: " << pose.position.transpose() << ", q: " << pose.rotation.w() << ", " << pose.rotation.x() << ", " << pose.rotation.y() << ", " << pose.rotation.z() << std::endl;


  // finally, transform the cloud:

  std::ofstream ofs;
  if (!ray::writePlyChunkStart(aligned_file, ofs))
    return false;

  // By maintaining these buffers below, we avoid almost all memory fragmentation  
  ray::RayPlyBuffer buffer;
  std::vector<Eigen::Vector3d> new_starts, new_ends;

  auto transform = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours)
  {
    for (size_t i = 0; i<ends.size(); i++)
    {
      new_starts.push_back(pose * starts[i]);
      new_ends.push_back(pose * ends[i]);
    }
    ray::writePlyChunk(ofs, buffer, new_starts, new_ends, times, colours);
    new_starts.clear();
    new_ends.clear();
  };

  if (!ray::readPly(cloud_name, true, transform, 0))
    return false;

  ray::writePlyChunkEnd(ofs);
  return true;  
}

} // namespace ray
