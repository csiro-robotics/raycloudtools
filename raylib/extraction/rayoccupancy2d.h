// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYOCCUPANCY2D_H
#define RAYLIB_RAYOCCUPANCY2D_H

#include "raylib/raylibconfig.h"
#include "raytrunks.h"
#include "../rayutils.h"

namespace ray
{
// grid structure as input to topology optimisation
struct Occupancy2D
{
  void init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width)
  {
    min_bound_ = min_bound;
    pixel_width_ = pixel_width;
    Eigen::Vector3d extent = max_bound - min_bound;
    dims_ = Eigen::Vector3d(std::ceil(extent[0] / pixel_width), std::ceil(extent[1] / pixel_width),
                            std::ceil(extent[2] / pixel_width))
              .cast<int>();
    std::cout << "min: " << min_bound.transpose() << ", ext: " << extent.transpose() << ", dims: " << dims_.transpose()
              << std::endl;

    pixels_.resize(dims_[0] * dims_[1]);
    memset(&pixels_[0], 0, sizeof(Pixel) * pixels_.size());
  }
  void save(const std::string &filename)
  {
    std::ofstream out(filename, std::ofstream::out);
    writePlainOldDataArray(out, pixels_);
    writePlainOldData(out, dims_);
    writePlainOldData(out, min_bound_);
    writePlainOldData(out, pixel_width_);
  }
  bool load(const std::string &filename)
  {
    std::ifstream input(filename, std::ifstream::in);
    if (!input.good())
      return false;
    readPlainOldDataArray(input, pixels_);
    readPlainOldData(input, dims_);
    readPlainOldData(input, min_bound_);
    readPlainOldData(input, pixel_width_);
    return true;
  }
  static constexpr int subpixels = 4;
  struct Pixel
  {
    uint16_t bits;
    inline double density() const { return 1.0 - ((double)bits / (double)(subpixels * subpixels)); }
  };
  Eigen::Vector3i dims_;
  Eigen::Vector3d min_bound_;
  double pixel_width_;
  std::vector<Pixel> pixels_;
  Pixel dummy_pixel_;
  inline Eigen::Vector3i pixelIndex(const Eigen::Vector3d &pos) const
  {
    return ((pos - min_bound_) / pixel_width_).cast<int>();
  }
  inline const Pixel &pixel(const Eigen::Vector3i &index) const
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline Pixel &pixel(const Eigen::Vector3i &index)
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline const Pixel &pixel(const Eigen::Vector3d &pos) const { return pixel(pixelIndex(pos)); }
  inline Pixel &pixel(const Eigen::Vector3d &pos) { return pixel(pixelIndex(pos)); }

  void fillDensities(const std::string &cloudname, const Eigen::ArrayXXd &lows, double clip_min, double clip_max)
  {
    ray::Cuboid bounds_;
    const double eps = 1e-9;
    bounds_.min_bound_ = min_bound_ + Eigen::Vector3d(eps, eps, eps);
    bounds_.max_bound_ = min_bound_ + dims_.cast<double>() * pixel_width_ - Eigen::Vector3d(eps, eps, eps);
    const double scale = (double)subpixels;

    auto addFreeSpace = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                            std::vector<double> &, std::vector<ray::RGBA> &) {
      for (size_t i = 0; i < ends.size(); ++i)
      {
        Eigen::Vector3d start = starts[i];
        Eigen::Vector3d end = ends[i];
        bounds_.clipRay(start, end);

        // now walk the pixels
        const Eigen::Vector3d dir = scale * (end - start);
        const Eigen::Vector3d source = scale * (start - min_bound_) / pixel_width_;
        const Eigen::Vector3d target = scale * (end - min_bound_) / pixel_width_;
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
          double ls[2] = { (round(p[0] + offsets[0]) - p[0]) / dir[0], (round(p[1] + offsets[1]) - p[1]) / dir[1] };
          int axis = (ls[0] < ls[1]) ? 0 : 1;
          inds[axis] += adds[axis];
          if (inds[axis] < 0 || inds[axis] >= subpixels * dims_[axis])
          {
            break;
          }
          double minL = ls[axis] * length;
          depth += minL + eps;
          p = source + dir * (depth / length);

          Eigen::Vector3i index = inds / subpixels;

          Eigen::Vector3d world_point = start + (end - start) * (depth / length);
          double height = world_point[2] - lows(index[0], index[1]);
          if (height > clip_min && height < clip_max)
          {
            Eigen::Vector3i rem = inds - subpixels * index;
            uint16_t bit = uint16_t(subpixels * rem[0] + rem[1]);
            pixel(index).bits |= uint16_t(1 << bit);
          }
        } while (depth <= maxDist);
      }
    };
    ray::Cloud::read(cloudname, addFreeSpace);
    //  draw("freespace_mess.png");

    auto removeOccupiedSpace = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                   std::vector<double> &, std::vector<ray::RGBA> &colours) {
      for (size_t i = 0; i < ends.size(); ++i)
      {
        if (colours[i].alpha == 0)
          continue;
//#define REMOVE_WHOLE_VOXEL
#if defined REMOVE_WHOLE_VOXEL
        const Eigen::Vector3d p = (ends[i] - min_bound_) / pixel_width_;
        const Eigen::Vector3i index = p.cast<int>();

        double height = ends[i][2] - lows(index[0], index[1]);
        if (height > 0.5 && height < clip_height)
          pixel(index).bits = 0;
#else
        const Eigen::Vector3d p = scale * (ends[i] - min_bound_) / pixel_width_;
        const Eigen::Vector3i inds = p.cast<int>();
        const Eigen::Vector3i index = inds / subpixels;
        const Eigen::Vector3i rem = inds - subpixels * index;
        uint16_t bit = uint16_t(subpixels * rem[0] + rem[1]);

        double height = ends[i][2] - lows(index[0], index[1]);
        if (height > clip_min && height < clip_max)
          pixel(index).bits &= (uint16_t) ~(uint16_t(1 << bit));
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
        if (vox.bits & ((uint16_t)1 << i))
          count++;
      vox.bits = count;
      bitcount += vox.bits;
    }
    //  draw("occupiedspace.png");

    std::cout << "average bit count: " << (double)bitcount / (double)pixels_.size() << std::endl;
  }

  void draw(const std::string &filename);
};


}  // namespace ray

#endif  // RAYLIB_RAYOCCUPANCY2D_H
