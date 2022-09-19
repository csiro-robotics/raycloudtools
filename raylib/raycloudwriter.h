// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUDWRITER_H
#define RAYLIB_RAYCLOUDWRITER_H

#include "raylib/raylibconfig.h"
#include "rayply.h"

namespace ray
{
/// This helper class is for writing a ray cloud to a file, one chunk at a time
/// These chunks can be any size, even 0
class RAYLIB_EXPORT CloudWriter
{
public:
  /// Open the file to write to
  bool begin(const std::string &file_name);

  /// write a set of rays to the file
  bool writeChunk(const class Cloud &chunk);

  /// write a set of rays to the file, direct arguments
  bool writeChunk(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times,
                  std::vector<RGBA> &colours)
  {
    return writeRayCloudChunk(ofs_, buffer_, starts, ends, times, colours, has_warned_);
  }

  /// finish writing, and adjust the vertex count at the start.
  void end();

  /// return the stored file name
  const std::string &fileName() { return file_name_; }

private:
  /// store the output file stream
  std::ofstream ofs_;
  /// store the file name, in order to provide a clear 'saved' message on end()
  std::string file_name_;
  /// ray buffer to avoid repeated reallocations
  RayPlyBuffer buffer_;
  /// whether a warning has been issued or not. This prevents multiple warnings.
  bool has_warned_;
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUDWRITER_H
