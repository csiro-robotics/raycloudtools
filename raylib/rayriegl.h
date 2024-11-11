// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tim Devereux

#ifndef RAYLIB_RAYRIEGL_H
#define RAYLIB_RAYRIEGL_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

#if RAYLIB_WITH_RIEGL
#include <riegl/scanlib.hpp>
#endif  // RAYLIB_WITH_RIEGL

namespace ray
{

/// Read an rxp file, into the fields passed by reference.
bool RAYLIB_EXPORT readRXP(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
                           std::vector<RGBA> &colours, double max_intensity,
                           std::vector<double> pose_transformation = std::vector<double>());

/// Chunk-based version of readRXP. This calls @c apply for every @c chunk_size points loaded
bool RAYLIB_EXPORT readRXP(const std::string &file_name,
                           std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<RGBA> &colours)>
                             apply,
                           size_t &num_bounded, double max_intensity, std::vector<double> pose_transformation,
                           size_t chunk_size = 1000000);

}  // namespace ray

#endif  // RAYLIB_RAYRIEGL_H