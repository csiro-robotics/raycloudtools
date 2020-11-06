// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLAZ_H
#define RAYLIB_RAYLAZ_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

namespace ray
{
/// Read a laz or las file, into the fields passed by reference. 
bool RAYLIB_EXPORT readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
                           std::vector<RGBA> &colours);

/// Chunk-based version of readLas. This calls @c apply for every @c chunk_size points loaded
bool RAYLIB_EXPORT readLas(const std::string &file_name,
     std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, 
     std::vector<double> &times, std::vector<RGBA> &colours)> apply, size_t &num_bounded, size_t chunk_size = 1000000);


/// Write to a laz or las file. The intensity is the only part that is extracted from the @c colours argument.
bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                            const std::vector<RGBA> &colours);
}

#endif  // RAYLIB_RAYLAZ_H
