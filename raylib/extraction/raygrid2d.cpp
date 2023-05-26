// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raygrid2d.h"

namespace ray
{
/// initialise for a given bounds and pixel width
void OccupancyGrid2D::init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width)
{
  min_bound_ = min_bound;
  pixel_width_ = pixel_width;
  const Eigen::Vector3d extent = max_bound - min_bound;
  dims_ = Eigen::Vector3d(std::ceil(extent[0] / pixel_width), std::ceil(extent[1] / pixel_width),
                          std::ceil(extent[2] / pixel_width))
            .cast<int>();
  std::cout << "min: " << min_bound.transpose() << ", ext: " << extent.transpose() << ", dims: " << dims_.transpose()
            << std::endl;

  pixels_.resize(dims_[0] * dims_[1]);
  memset(&pixels_[0], 0, sizeof(Pixel) * pixels_.size());
}

/// save the grid
void OccupancyGrid2D::save(const std::string &filename)
{
  std::ofstream out(filename, std::ofstream::out);
  writePlainOldDataArray(out, pixels_);
  writePlainOldData(out, dims_);
  writePlainOldData(out, min_bound_);
  writePlainOldData(out, pixel_width_);
}

/// load the grid
bool OccupancyGrid2D::load(const std::string &filename)
{
  std::ifstream input(filename, std::ifstream::in);
  if (!input.good())
  {
    return false;
  }
  readPlainOldDataArray(input, pixels_);
  readPlainOldData(input, dims_);
  readPlainOldData(input, min_bound_);
  readPlainOldData(input, pixel_width_);
  return true;
}

/// fill in the occupancy data (pixels_) based on the rays in the cloud @c cloudname,
/// within a height window @c clip_min to @c clip_max
void OccupancyGrid2D::fillDensities(const std::string &cloudname, const Eigen::ArrayXXd &lows, double clip_min,
                                    double clip_max)
{
  ray::Cuboid bounds_;
  const double eps = 1e-9;
  bounds_.min_bound_ = min_bound_ + Eigen::Vector3d(eps, eps, eps);
  bounds_.max_bound_ = min_bound_ + dims_.cast<double>() * pixel_width_ - Eigen::Vector3d(eps, eps, eps);
  const double scale = static_cast<double>(GRID2D_SUBPIXELS);

  // filling in the free space per chunk of ray cloud
  auto addFreeSpace = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                          std::vector<double> &, std::vector<ray::RGBA> &) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end)) // clip the ray within the bounds
      {
        continue;
      }  

      // now walk the pixels
      const Eigen::Vector3d dir = scale * (end - start);
      const Eigen::Vector3d source = scale * (start - min_bound_) / pixel_width_;
      const Eigen::Vector3d target = scale * (end - min_bound_) / pixel_width_;
      const double length = dir.norm();
      const double eps = 1e-9;  // to stay away from edge cases
      // remove 2 GRID2D_SUBPIXELS to give a small buffer around the object
      const double maxDist = (target - source).norm() - 2.0;

      // cached values to speed up the loop below
      Eigen::Vector3i adds;
      Eigen::Vector3d offsets;
      for (int k = 0; k < 3; ++k)
      {
        if (dir[k] > 0.0)
        {
          adds[k] = 1;
          offsets[k] = 0.5;
        }
        else
        {
          adds[k] = -1;
          offsets[k] = -0.5;
        }
      }

      Eigen::Vector3d p = source;  // our moving variable as we walk over the grid
      Eigen::Vector3i inds = p.cast<int>();
      double depth = 0;
      // walk over the grid, one pixel at a time.
      do
      {
        // deltas in each axis
        const double ls[2] = { (round(p[0] + offsets[0]) - p[0]) / dir[0], (round(p[1] + offsets[1]) - p[1]) / dir[1] };
        // shift in the axis with smallest delta
        int axis = (ls[0] < ls[1]) ? 0 : 1;
        // update the index to the new cell
        inds[axis] += adds[axis];
        if (inds[axis] < 0 || inds[axis] >= GRID2D_SUBPIXELS * dims_[axis])
        {
          break;
        }
        // minimum length of line segment within cell
        const double minL = ls[axis] * length;
        depth += minL + eps;
        // update position p
        p = source + dir * (depth / length);

        // get the index of the pixel
        Eigen::Vector3i index = inds / GRID2D_SUBPIXELS;

        // find the world space location
        Eigen::Vector3d world_point = start + (end - start) * (depth / length);
        // get the height above ground at this location
        const double height = world_point[2] - lows(index[0], index[1]);
        if (height > clip_min && height < clip_max)  // only update occupancy within height window
        {
          // some bit trickery to fill in part of the 4x4 grid per pixel
          const Eigen::Vector3i rem = inds - GRID2D_SUBPIXELS * index;
          const uint16_t bit = uint16_t(GRID2D_SUBPIXELS * rem[0] + rem[1]);
          pixel(index).bits |= uint16_t(1 << bit);
        }
      } while (depth <= maxDist);
    }
  };
  ray::Cloud::read(cloudname, addFreeSpace);

  // wherever these is an end point, we want to remove it as free space
  auto removeOccupiedSpace = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                 std::vector<double> &, std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      if (colours[i].alpha == 0)
      {
        continue;
      }
//#define REMOVE_WHOLE_VOXEL
#if defined REMOVE_WHOLE_VOXEL
      const Eigen::Vector3d p = (ends[i] - min_bound_) / pixel_width_;
      const Eigen::Vector3i index = p.cast<int>();

      double height = ends[i][2] - lows(index[0], index[1]);
      if (height > 0.5 && height < clip_height)
        pixel(index).bits = 0;
#else
      // find the subpixel that this point is in
      const Eigen::Vector3d p = scale * (ends[i] - min_bound_) / pixel_width_;
      const Eigen::Vector3i inds = p.cast<int>();
      const Eigen::Vector3i index = inds / GRID2D_SUBPIXELS;
      const Eigen::Vector3i rem = inds - GRID2D_SUBPIXELS * index;
      uint16_t bit = uint16_t(GRID2D_SUBPIXELS * rem[0] + rem[1]);

      double height = ends[i][2] - lows(index[0], index[1]);
      if (height > clip_min && height < clip_max)               // if within the height window
        pixel(index).bits &= (uint16_t) ~(uint16_t(1 << bit));  // then remove it
#endif
    }
  };
  ray::Cloud::read(cloudname, removeOccupiedSpace);

  // convert the bit fields into subpixel counts
  unsigned long bitcount = 0;
  for (auto &vox : pixels_)
  {
    uint16_t count = 0;
    for (unsigned long i = 0; i < 16; i++)
    {
      if (vox.bits & ((uint16_t)1 << i))
      {
        count++;
      }
    }
    vox.bits = count;
    bitcount += vox.bits;
  }

  std::cout << "average bit count: " << static_cast<double>(bitcount) / static_cast<double>(pixels_.size())
            << std::endl;
}

void RayIndexGrid2D::init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width)
{
  min_bound_ = min_bound;
  pixel_width_ = pixel_width;
  const Eigen::Vector3d extent = max_bound - min_bound;
  dims_ = Eigen::Vector3d(std::ceil(extent[0] / pixel_width), std::ceil(extent[1] / pixel_width),
                          std::ceil(extent[2] / pixel_width))
            .cast<int>();
  std::cout << "min: " << min_bound.transpose() << ", ext: " << extent.transpose() << ", dims: " << dims_.transpose()
            << std::endl;

  pixels_.resize(dims_[0] * dims_[1]);
  memset(&pixels_[0], 0, sizeof(Pixel) * pixels_.size());
}

void RayIndexGrid2D::fillRays(const Cloud &cloud)
{
  ray::Cuboid bounds_;
  const double eps = 1e-9;
  bounds_.min_bound_ = min_bound_ + Eigen::Vector3d(eps, eps, eps);
  bounds_.max_bound_ = min_bound_ + dims_.cast<double>() * pixel_width_ - Eigen::Vector3d(eps, eps, eps);

  for (size_t i = 0; i < cloud.ends.size(); ++i)
  {
    Eigen::Vector3d start = cloud.starts[i];
    Eigen::Vector3d end = cloud.ends[i];
    if (!bounds_.clipRay(start, end))
    {
      continue;
    }

    // now walk the pixels
    const Eigen::Vector3d dir = end - start;
    const Eigen::Vector3d source = (start - min_bound_) / pixel_width_;
    const Eigen::Vector3d target = (end - min_bound_) / pixel_width_;
    const double length = dir.norm();
    const double eps = 1e-9;  // to stay away from edge cases
    const double maxDist =
      (target - source).norm() - 2.0;  // remove 2 subpixels to give a small buffer around the object

    // cached values to speed up the loop below
    Eigen::Vector3i adds;
    Eigen::Vector3d offsets;
    for (int k = 0; k < 3; ++k)
    {
      if (dir[k] > 0.0)
      {
        adds[k] = 1;
        offsets[k] = 0.5;
      }
      else
      {
        adds[k] = -1;
        offsets[k] = -0.5;
      }
    }

    Eigen::Vector3d p = source;  // our moving variable as we walk over the grid
    Eigen::Vector3i inds = p.cast<int>();
    double depth = 0;
    // walk over the grid, one pixel at a time.
    do
    {
      const double ls[2] = { (round(p[0] + offsets[0]) - p[0]) / dir[0], (round(p[1] + offsets[1]) - p[1]) / dir[1] };
      const int axis = (ls[0] < ls[1]) ? 0 : 1;
      inds[axis] += adds[axis];
      if (inds[axis] < 0 || inds[axis] >= dims_[axis])
      {
        break;
      }
      const double minL = ls[axis] * length;
      depth += minL + eps;
      p = source + dir * (depth / length);
      Pixel &pix = pixel(inds);
      if (pix.filled)
      {
        pix.ray_ids.push_back(static_cast<int>(i));
      }
    } while (depth <= maxDist);
  }
}
}  // namespace ray
