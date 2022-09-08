// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLAZ_H
#define RAYLIB_RAYLAZ_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

#if RAYLIB_WITH_LAS
#include <liblas/reader.hpp>
#endif  // RAYLIB_WITH_LAS


namespace ray
{
/// Read a laz or las file, into the fields passed by reference.
bool RAYLIB_EXPORT readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
                           std::vector<RGBA> &colours, double max_intensity,
                           Eigen::Vector3d *offset_to_remove = nullptr);

/// Chunk-based version of readLas. This calls @c apply for every @c chunk_size points loaded
bool RAYLIB_EXPORT readLas(const std::string &file_name,
                           std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<RGBA> &colours)>
                             apply,
                           size_t &num_bounded, double max_intensity, Eigen::Vector3d *offset_to_remove,
                           size_t chunk_size = 1000000);


/// Write to a laz or las file. The intensity is the only part that is extracted from the @c colours argument.
bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points,
                            const std::vector<double> &times, const std::vector<RGBA> &colours);

/// Class for chunked writing of las/laz files.
class RAYLIB_EXPORT LasWriter
{
public:
  /// construct the class with a file name, which is stored
  LasWriter(const std::string &file_name);
  /// the destructor
  ~LasWriter();
  /// write a chunk of points to the file, described by the vector arguments
  bool writeChunk(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                  const std::vector<RGBA> &colours);

private:
  const std::string &file_name_;
  std::ofstream out_;
#if RAYLIB_WITH_LAS
  liblas::Header header_;
  liblas::Writer *writer_;
#endif  // RAYLIB_WITH_LAS
};
}  // namespace ray

#endif  // RAYLIB_RAYLAZ_H
