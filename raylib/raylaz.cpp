// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylaz.h"

#include "rayunused.h"

#if RAYLIB_WITH_LAS
#include <liblas/reader.hpp>
#include <liblas/factory.hpp>
#include <liblas/point.hpp>
#endif  // RAYLIB_WITH_LAS

namespace ray
{
bool readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times, std::vector<RGBA> &colours,
                  int decimation)
{
#if RAYLIB_WITH_LAS
  std::vector<uint8_t> intensities;
  std::cout << "readLas: filename: " << file_name << std::endl;

  std::ifstream ifs;
  ifs.open(file_name.c_str(), std::ios::in | std::ios::binary);

  if (!ifs.is_open())
  {
    std::cerr << "readLas: failed to open stream" << std::endl;
    return false;
  }

  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  liblas::Header header = reader.GetHeader();

  unsigned int size = header.GetPointRecordsCount();

  int num_intensities = 0;
  for (unsigned int i = 0; i < size; i++)
  {
    reader.ReadNextPoint();
    liblas::Point point = reader.GetPoint();
    // we're downsampling to 1/decimation of the orginal file.
    if ((i % decimation) == 0)
    {
      Eigen::Vector3d position;
      position[0] = point.GetX();
      position[1] = point.GetY();
      position[2] = point.GetZ();
      positions.push_back(position);
      times.push_back(point.GetTime());
      uint8_t intensity = (uint8_t) std::min(point.GetIntensity(), (uint16_t)255);
      if (intensity > 0)
        num_intensities++;
      intensities.push_back(intensity);
    }
  }
  if (num_intensities == 0)
    for (auto &i : intensities) i = 1.0;
  colourByTime(times, colours);
  for (int i = 0; i < (int)colours.size(); i++)  // add intensity into alhpa channel
    colours[i].alpha = intensities[i];

  std::cout << "loaded " << file_name << " with " << positions.size() << " points" << std::endl;
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(positions);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  RAYLIB_UNUSED(decimation);
  std::cerr << "readLas: cannot read file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                            const std::vector<RGBA> &colours)
{
#if RAYLIB_WITH_LAS
  std::cout << "saving LAZ file" << std::endl;
 
  liblas::Header header;
  header.SetDataFormatId(liblas::ePointFormat1); // Time only

  if (file_name.find(".laz") != std::string::npos)
    header.SetCompressed(true);
 
  std::cout << "Saving points to " << file_name << std::endl;

  std::ofstream ofs;
  ofs.open(file_name.c_str(), std::ios::out | std::ios::binary);
  if (!ofs.is_open())
    return false;

  const double scale = 1e-4;
  header.SetScale(scale, scale, scale);

  liblas::Writer writer(ofs, header);

  liblas::Point point(&header);
  point.SetHeader(&header);//TODO HACK Version 1.7.0 does not correctly resize the data. Commit 6e8657336ba445fcec3c9e70c2ebcd2e25af40b9 (1.8.0 3 July fixes it)
  for (unsigned int i = 0; i < points.size(); i++)
  {
    point.SetCoordinates(points[i][0], points[i][1], points[i][2]);
    point.SetIntensity(colours[i].alpha);
    if (!times.empty())
      point.SetTime(times[i]);
    writer.WritePoint(point);
  }
  return true;
#else // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif // RAYLIB_WITH_LAS
}

} // ray