// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayaxisalign.h"
#include "raycloud.h"
#include "raycloudwriter.h"
#include "rayply.h"
#include "rayunused.h"
#include "rayutils.h"

namespace ray
{
namespace
{
/// These constants are the angular and amplitude (distance from centre) resolutions for the radon transform.
/// Too small and the alignment is imprecise, too large and it fails to align on noisy planes (such as a crop row).
/// We may ultimately make these a function of the cloud's size or density, but for this iteration of the axis alignment
/// function, they are fixed.
const int ang_res = 256;  // ang_res must be divisible by 2
const int amp_res = 256;
}  // namespace

/// the output file @c out_file is the input file @c in_file, transformed by @c pose
/// returns whether the cloud was successfully modified
bool transformAndSaveCloud(const std::string &in_file, const std::string &out_file, const Pose &pose)
{
  CloudWriter writer;
  if (!writer.begin(out_file))
    return false;
  Cloud chunk;

  auto transform = [&chunk, &writer, &pose](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                            std::vector<double> &times, std::vector<RGBA> &colours) {
    chunk.clear();
    chunk.times = times;
    chunk.colours = colours;
    for (size_t i = 0; i < ends.size(); i++)
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

// Apply the radon transform to the 2D array representing density of end points (and their centroid), and convert the
// peak into a 2D pose For convenience we use the 3D Euclidean transformation class (Pose), and set the vertical
// translation to the mid point.
Pose estimate2DPose(const Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> &position_accumulator,
                    const Cloud::Info &info)
{
  const Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  const Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;
  const Eigen::Vector3d mid_bound = (min_bound + max_bound) / 2.0;

  // radius used for maximum amplitude (distance from mid_bound)
  const double radius = 0.5 * std::sqrt(ray::sqr(max_bound[0] - min_bound[0]) + ray::sqr(max_bound[1] - min_bound[1]));
  const double eps = 0.0001;  // avoids the most distant point exceeding the array bounds

  Eigen::ArrayXXd weights(ang_res, amp_res);  // this is the output of the radon transform
  weights.fill(0);
  // Radon transform
  for (int ii = 0; ii < amp_res; ii++)
  {
    for (int jj = 0; jj < amp_res; jj++)
    {
      const Eigen::Vector3d &accumulator = position_accumulator(ii, jj);
      const double weight =
        accumulator[2];  // the weight here is the number of end points under this pixel (array cell)
      if (weight == 0.0)
        continue;
      const Eigen::Vector2d centroid(accumulator[0] / weight, accumulator[1] / weight);
      const double angle = atan2(centroid[0], centroid[1]);
      const double amplitude = std::sqrt(centroid[0] * centroid[0] + centroid[1] * centroid[1]) / radius;

      // now draw the sine wave for this point.
      for (int i = 0; i < ang_res; i++)
      {
        const double ang = kPi * static_cast<double>(i) / static_cast<double>(ang_res);
        const double height = amplitude * std::sin(ang + angle);
        const double y = (static_cast<double>(amp_res) - 1.0 - eps) * (0.5 + 0.5 * height);  // rescale the sine wave

        // linear blend of the weight onto the two nearest neighbour pixels
        const int j = static_cast<int>(y);
        const double blend = y - static_cast<double>(j);
        weights(i, j) += (1.0 - blend) * weight;
        weights(i, j + 1) += blend * weight;
      }
    }
  }
  // now find greatest weight cell:
  int max_i = 0, max_j = 0;
  weights.maxCoeff(&max_i, &max_j);

  // find the sub-pixel peak, using its two neighbours in angle ...
  double peak_angle = static_cast<double>(max_i) + peak(weights((max_i + ang_res - 1) % ang_res, max_j),
                                                        weights(max_i, max_j), weights((max_i + 1) % ang_res, max_j));
  peak_angle *= kPi / static_cast<double>(ang_res);
  std::cout << "principle angle: " << peak_angle << ", direction vector: " << std::cos(peak_angle) << ", "
            << std::sin(peak_angle) << std::endl;
  // ... and its two neighbours in the amplitude axis
  double peak_amp = static_cast<double>(max_j) + peak(weights(max_i, std::max(0, max_j - 1)), weights(max_i, max_j),
                                                      weights(max_i, std::min(max_j + 1, amp_res - 1)));
  peak_amp = radius * ((2.0 * peak_amp / static_cast<double>(amp_res)) - 1.0);
  std::cout << "distance along direction: " << peak_amp << std::endl;

  // vector representing the strongest plane
  Eigen::Vector2d line_vector = peak_amp * Eigen::Vector2d(std::cos(peak_angle), std::sin(peak_angle));

  // now find the orthogonal best edge, which is the greatest weight along the orthogonal angle
  const int orth_i = (max_i + ang_res / 2) % ang_res;
  int max_orth_j = 0;
  weights.row(orth_i).maxCoeff(&max_orth_j);

  // now find the sub-pixel peak amplitude in this orthogonal direction:
  double ortho_amp =
    static_cast<double>(max_orth_j) + peak(weights(orth_i, std::max(0, max_orth_j - 1)), weights(orth_i, max_orth_j),
                                           weights(orth_i, std::min(max_orth_j + 1, amp_res - 1)));
  ortho_amp = radius * ((2.0 * ortho_amp / static_cast<double>(amp_res)) - 1.0);
  if (orth_i < max_i)        // the %res earlier puts it in antiphase (since we only have 180 degrees per weights map)
    ortho_amp = -ortho_amp;  // so negate the amplitude here

  const Eigen::Vector2d orthogonal_line_vector =
    ortho_amp * Eigen::Vector2d(std::cos(peak_angle + kPi / 2.0), std::sin(peak_angle + kPi / 2.0));
  std::cout << "primary amplitude: " << peak_amp << ", orthogonal amplitude: " << ortho_amp << std::endl;

  // we can now find the 2D pose of best fit
  const Eigen::Vector2d centre =
    line_vector + orthogonal_line_vector;  // find the intersection of the orthogonal strongest planes
  const Pose to_mid(-mid_bound - Eigen::Vector3d(centre[0], centre[1], 0), Eigen::Quaterniond::Identity());
  const Pose rotation(Eigen::Vector3d::Zero(),
                      Eigen::Quaterniond(Eigen::AngleAxisd(-peak_angle, Eigen::Vector3d(0, 0, 1))));
  Pose pose = rotation * to_mid;

  // now rotate the cloud into the positive x/y axis, depending on which has the largest difference
  const Eigen::Vector3d mid_point = pose * info.centroid;
  const Eigen::Quaterniond yaw180(0, 0, 0, 1);
  if (std::abs(mid_point[0]) > std::abs(mid_point[1]))
  {
    if (mid_point[0] < 0.0)
      pose = Pose(Eigen::Vector3d::Zero(), yaw180) * pose;
  }
  else
  {
    if (mid_point[1] < 0.0)
      pose = Pose(Eigen::Vector3d::Zero(), yaw180) * pose;
  }
  return pose;
}

// A radon transform is used as follows:
// 1. we quantise the cloud into a 2D grid of centroids, weighted by number of end points within the cell
// 2. for each cell, render a sine wave to weights image
// 3. find the maximum weight in the image (the new y axis) and interpolate the angle and position using its neighbours
// 4. find the highest orthogonal weight (the new x axis) and interpolate its position
// 5. quantise density vertically into an array to get the strongest ground height signal, interpolating the max value
// 6. transform the cloud according to these axes, and save it to file
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
  const Eigen::Vector3d mid_bound = (min_bound + max_bound) / 2.0;

  // fill in the arrays
  double eps = 0.0001;  // to stop edge cases exceeding the array bounds

  // 1. Convert the cloud into a weighted centroid field. I'm using element [2] for the weight.
  Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> position_accumulator(amp_res, amp_res);
  position_accumulator.fill(Eigen::Vector3d::Zero());
  double step_x = (static_cast<double>(amp_res) - 1.0 - eps) / (max_bound[0] - min_bound[0]);
  double step_y = (static_cast<double>(amp_res) - 1.0 - eps) / (max_bound[1] - min_bound[1]);
  double step_z = (static_cast<double>(amp_res) - 1.0 - eps) / (max_bound[2] - min_bound[2]);
  // and set up the vertical arrays too
  Eigen::ArrayXd vertical_weights(
    amp_res);  // TODO: height is usually much less than width... more constant voxel size?
  vertical_weights.fill(0);

  // fill in the accumulator: (sum of positions, #end points) within each cell
  // also fill in the vertical density (number of end points) array
  auto fill_arrays = [&min_bound, &mid_bound, &position_accumulator, &step_x, &step_y, &step_z, &vertical_weights](
                       std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    RAYLIB_UNUSED(starts);
    RAYLIB_UNUSED(times);
    for (size_t e = 0; e < ends.size(); e++)
    {
      if (colours[e].alpha == 0)  // unbounded
        continue;
      Eigen::Vector3d index = ends[e] - min_bound;  // find the cell index for this ray end-point
      index[0] *= step_x;
      index[1] *= step_y;
      index[2] *= step_z;
      Eigen::Vector3d pos = ends[e] - mid_bound;
      pos[2] = 1.0;  // this element counts the number of end points within the cell
      position_accumulator((int)index[0], (int)index[1]) += pos;

      // Distribute the weighting linearly between two nearest cells in the 1D array
      const int k = (int)index[2];
      const double blend = index[2] - static_cast<double>(k);
      vertical_weights[k] += (1.0 - blend);
      vertical_weights[k + 1] += blend;
    }
  };
  if (!Cloud::read(cloud_name, fill_arrays))
    return false;

  // 2,3,4. Apply the radon transform in 2D
  Pose pose = estimate2DPose(position_accumulator, info);

  // 5. get the vertical displacement
  int max_k = 0;  // k is the vertical cell index, as in (i,j,k)
  vertical_weights.maxCoeff(&max_k);

  double peak_k = static_cast<double>(max_k) + peak(vertical_weights[std::max(0, max_k - 1)], vertical_weights[max_k],
                                                    vertical_weights[std::min(max_k + 1, amp_res - 1)]);
  pose.position[2] = -(peak_k / step_z) - min_bound[2];

  std::cout << "pose: " << pose.position.transpose() << ", q: " << pose.rotation.w() << ", " << pose.rotation.x()
            << ", " << pose.rotation.y() << ", " << pose.rotation.z() << std::endl;

  // 6. transform the cloud and save to file
  return transformAndSaveCloud(cloud_name, aligned_file, pose);
}

}  // namespace ray
