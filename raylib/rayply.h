// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPLY_H
#define RAYLIB_RAYPLY_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

namespace ray
{
#if RAYLIB_DOUBLE_RAYS
using PointPlyEntry = Eigen::Matrix<float, 9, 1>;
using RayPlyEntry = Eigen::Matrix<float, 12, 1>;  // structure of raycloud cloud rays, written to ply file
#else
using PointPlyEntry = Eigen::Matrix<float, 6, 1>;
using RayPlyEntry = Eigen::Matrix<float, 9, 1>;  // structure of raycloud cloud rays, written to ply file
#endif
using PointPlyBuffer = std::vector<PointPlyEntry>;
using RayPlyBuffer = std::vector<RayPlyEntry>;  // buffer for storing a list of rays to be written

/// read in a .ply file into the fields given by reference
/// Note that @c max_intensity is only used when reading in a point cloud. Intensities are already stored in the
/// colour alpha channel in ray clouds.
bool RAYLIB_EXPORT readPly(const std::string &file_name, std::vector<Eigen::Vector3d> &starts,
                           std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours,
                           bool is_ray_cloud, double max_intensity = 0);
/// read in a .ply file that represents a triangular mesh, into the @c Mesh structure
bool RAYLIB_EXPORT readPlyMesh(const std::string &file, class Mesh &mesh);

/// write a .ply file representing a triangular mesh
bool RAYLIB_EXPORT writePlyMesh(const std::string &file_name, const class Mesh &mesh, bool flip_normals = false);

/// ready in a ray cloud or point cloud .ply file, and call the @c apply function one chunk at a time,
/// @c chunk_size is the number of rays to read at one time. This method can be used on large clouds where
/// the full set of rays is not required to be in memory at one time.
/// @c times_optional flag allows clouds to be read with no time stamps
bool RAYLIB_EXPORT readPly(const std::string &file_name, bool is_ray_cloud,
                           std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<RGBA> &colours)>
                             apply, 
                           double max_intensity, bool times_optional = false, size_t chunk_size = 1000000);


/// write a .ply file representing a point cloud
bool RAYLIB_EXPORT writePlyPointCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &points,
                                      const std::vector<double> &times, const std::vector<RGBA> &colours);

/// write a .ply file representing a ray cloud
bool RAYLIB_EXPORT writePlyRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                                    const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                    const std::vector<RGBA> &colours);

/// Chunked version of writePlyRayCloud
bool RAYLIB_EXPORT writeRayCloudChunkStart(const std::string &file_name, std::ofstream &out);
bool RAYLIB_EXPORT writeRayCloudChunk(std::ofstream &out, RayPlyBuffer &vertices,
                                      const std::vector<Eigen::Vector3d> &starts,
                                      const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                      const std::vector<RGBA> &colours, bool &has_warned);
unsigned long RAYLIB_EXPORT writeRayCloudChunkEnd(std::ofstream &out);

/// Chunked version of writePlyPointCloud
bool RAYLIB_EXPORT writePointCloudChunkStart(const std::string &file_name, std::ofstream &out);
bool RAYLIB_EXPORT writePointCloudChunk(std::ofstream &out, PointPlyBuffer &vertices,
                                        const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                                        const std::vector<RGBA> &colours, bool &has_warned);
void RAYLIB_EXPORT writePointCloudChunkEnd(std::ofstream &out);

/// Simple function for converting a ray cloud according to the per-ray function @c apply
bool convertCloud(const std::string &in_name, const std::string &out_name,
                  std::function<void(Eigen::Vector3d &start, Eigen::Vector3d &ends, double &time, RGBA &colour)> apply);
}  // namespace ray

#endif  // RAYLIB_RAYPLY_H
