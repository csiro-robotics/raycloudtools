// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloudwriter.h"
#include "raycloud.h"

namespace ray
{

bool CloudWriter::begin(const std::string &file_name)
{
  if (file_name.empty())
  {
    std::cerr << "Error: cloud writer begin called with empty file name" << std::endl;
    return false;
  }
  file_name_ = file_name;
  if (!writePlyChunkStart(file_name_, ofs_))
  {
    std::cerr << "cannot write to file: " << file_name_ << std::endl;
    return false;    
  }
  return true;
}

void CloudWriter::end()
{
  if (file_name_.empty()) // no effect if begin has not been called
    return;
  const unsigned long num_rays = ray::writePlyChunkEnd(ofs_);
  std::cout << num_rays << " rays saved to " << file_name_ << std::endl;
  ofs_.close();
}

bool CloudWriter::writeChunk(const Cloud &chunk)
{
  return writePlyChunk(ofs_, buffer_, chunk.starts, chunk.ends, chunk.times, chunk.colours);
}


} // namespace ray
