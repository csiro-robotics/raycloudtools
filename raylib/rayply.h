// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPLY_H
#define RAYLIB_RAYPLY_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"
#include "raymesh.h"

namespace ray
{
/// read in a .ply file into the fields given by reference
bool RAYLIB_EXPORT readPly(const std::string &file_name, std::vector<Eigen::Vector3d> &starts,
                           std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours, 
                           bool is_ray_cloud);
/// read in a .ply file that represents a triangular mesh, into the @c Mesh structure
bool RAYLIB_EXPORT readPlyMesh(const std::string &file, Mesh &mesh);

/// write a .ply file representing a triangular mesh
void RAYLIB_EXPORT writePlyMesh(const std::string &file_name_rawaw, const Mesh &mesh, bool flip_normals = false);


bool readPly(const std::string &file_name, bool is_ray_cloud, 
     std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, 
     std::vector<double> &times, std::vector<RGBA> &colours)> apply, size_t chunk_size = 1000000);


/// write a .ply file representing a point cloud
void RAYLIB_EXPORT writePly(const std::string &file_name, const std::vector<Eigen::Vector3d> &points, 
                            const std::vector<double> &times, const std::vector<RGBA> &colours);

/// write a .ply file representing a ray cloud
void RAYLIB_EXPORT writePly(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                            const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                            const std::vector<RGBA> &colours);

/// Chunked version of writePly for ray clouds
bool RAYLIB_EXPORT writePlyChunkStart(const std::string &file_name, std::ofstream &out);
bool RAYLIB_EXPORT writePlyChunk(std::ofstream &out, std::vector<Eigen::Matrix<float, 9, 1>> &vertices, const std::vector<Eigen::Vector3d> &starts,
     const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times, const std::vector<RGBA> &colours);
void RAYLIB_EXPORT writePlyChunkEnd(std::ofstream &out);
}  // namespace ray

#endif  // RAYLIB_RAYPLY_H
