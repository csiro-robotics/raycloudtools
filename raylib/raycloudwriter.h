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

  /// finish writing, and adjust the vertex count at the start.
  void end();
private:
  RayPlyBuffer buffer_;
  std::ofstream ofs_;
  std::string file_name_;
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUDWRITER_H
