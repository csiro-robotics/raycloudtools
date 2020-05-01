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
bool RAYLIB_EXPORT readPly(const std::string &file_namee, std::vector<Eigen::Vector3d> &starts,
                           std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours);
bool RAYLIB_EXPORT readPlyMesh(const std::string &file, Mesh &mesh);

void RAYLIB_EXPORT writePly(const std::string &file_namee, const std::vector<Eigen::Vector3d> &starts,
                            const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                            const std::vector<RGBA> &colours);
void RAYLIB_EXPORT writePlyMesh(const std::string &file_name_rawaw, const Mesh &mesh, bool flip_normalss = false);
}  // namespace ray

#endif  // RAYLIB_RAYPLY_H
