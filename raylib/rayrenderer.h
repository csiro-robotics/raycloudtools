// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYRENDERER_H
#define RAYLIB_RAYRENDERER_H

#include "raylib/raylibconfig.h"
#include "raycuboid.h"
#include "rayutils.h"
#include "raypose.h"

namespace ray
{
/// Supported view directions on cloud data
enum class RAYLIB_EXPORT ViewDirection
{
  Top, 
  Left, 
  Right, 
  Front,
  Back
};

/// Supported render styles on cloud data
enum class RAYLIB_EXPORT RenderStyle
{
  Ends, 
  Mean, 
  Sum, 
  Starts, 
  Rays, 
  Density, 
  Density_rgb
};

/// Render a ray cloud according to the supplied parameters
bool RAYLIB_EXPORT renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction, 
                               RenderStyle style, double pix_width, const std::string &image_file, const std::string &projection_file);

/// This is used for estimating the per-voxel density of a ray cloud
/// Density represents the surface area per volume, assuming an unbiased distribution of surface angles
/// It is most effective as a measure of leaf area per volume on vegetation, and is described in:
/// Lowe, Thomas, et al. "Canopy Density Estimation in Perennial Horticulture Crops Using 3D Spinning LiDAR SLAM." 
/// arXiv preprint arXiv:2007.15652 (2020).
struct RAYLIB_EXPORT DensityGrid
{
  static const int min_voxel_hits = 2;
  static constexpr double spherical_distribution_scale = 2.0; // average area scale due to a spherical uniform distribution of leave angles relative to the rays
  
  DensityGrid(const Cuboid &grid_bounds, double vox_width, const Eigen::Vector3i &dims) : 
    bounds_(grid_bounds), voxel_width_(vox_width), voxel_dims_(dims) 
  {
    voxels_.resize(dims[0]*dims[1]*dims[2]);
  }
  
  /// This specific voxel class represents a density
  class Voxel
  {
  public:
    Voxel(){ num_hits_ = num_rays_ = path_length_ = 0.0; }
    /// return the density that the voxel represents
    inline double density() const;
    /// the densities can be summed element-wise
    inline void operator +=(const Voxel &other);
    /// the densities can be multiplied by a scalar, element-wise
    inline Voxel operator *(float scale) const;
    /// Add a ray which enters the voxel and hits within it. @c length is the ray length within the voxel
    inline void addHitRay(float length);
    /// Add a ray which enters and exits the voxel. @c length is the ray path length within the voxel
    inline void addMissRay(float length);
    inline const float &numHits() const { return num_hits_; }
    inline const float &numRays() const { return num_rays_; }
    inline const float &pathLength() const { return path_length_; }
    
  private:
    float num_hits_;
    float num_rays_;
    float path_length_;
  };

  /// This streams in a ray cloud file, and fills in the voxel density information
  void calculateDensities(const std::string &file_name);
  /// To void low-ray-count voxels giving unstable density estimates, we fuse with neighbour information
  /// up to a specified minimum number of rays. Specified in DENSITY_MIN_RAYS 
  void addNeighbourPriors();
  /// Note, for performance, this index function does not check that the specified indices are in valid bounds. 
  /// It is up to the calling function to assure this condition
  inline int getIndex(const Eigen::Vector3i &inds) const;
  /// Return the vector of density voxels
  inline const std::vector<Voxel> &voxels() const { return voxels_; }

private:
  Cuboid bounds_;
  std::vector<Voxel> voxels_;
  double voxel_width_;
  Eigen::Vector3i voxel_dims_;
};

// inline functions
inline double DensityGrid::Voxel::density() const 
{ 
  if (num_rays_ <= min_voxel_hits)
    return 0.0;
  return spherical_distribution_scale * (num_rays_-1.0) * num_hits_ / (1e-10 + num_rays_*path_length_); 
} 
void DensityGrid::Voxel::operator +=(const DensityGrid::Voxel &other)
{
  num_hits_ += other.num_hits_;
  num_rays_ += other.num_rays_; 
  path_length_ += other.path_length_; 
}
DensityGrid::Voxel DensityGrid::Voxel::operator *(float scale) const
{
  DensityGrid::Voxel voxel;
  voxel.num_hits_ = num_hits_ * scale;
  voxel.num_rays_ = num_rays_ * scale;
  voxel.path_length_ = path_length_ * scale;
  return voxel;
}
void DensityGrid::Voxel::addHitRay(float length)
{
  path_length_ += length;
  num_hits_++;
  num_rays_++;
}
void DensityGrid::Voxel::addMissRay(float length)
{
  path_length_ += length;
  num_rays_++;
}
int DensityGrid::getIndex(const Eigen::Vector3i &inds) const
{
  return inds[0] + inds[1]*voxel_dims_[0] + inds[2] * voxel_dims_[0]*voxel_dims_[1];
}
}
#endif // RAYLIB_RAYRENDERER_H
