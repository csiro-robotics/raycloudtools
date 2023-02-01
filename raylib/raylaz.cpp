// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylaz.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "rayunused.h"

#if RAYLIB_WITH_LAS
#include <liblas/factory.hpp>
#include <liblas/point.hpp>
#endif  // RAYLIB_WITH_LAS

namespace ray
{
bool readLas(const std::string &file_name,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply,
             size_t &num_bounded, double max_intensity, Eigen::Vector3d *offset_to_remove, size_t chunk_size)
{
#if RAYLIB_WITH_LAS
  std::cout << "readLas: filename: " << file_name << std::endl;

  std::ifstream ifs;
  ifs.open(file_name.c_str(), std::ios::in | std::ios::binary);

  if (ifs.fail())
  {
    std::cerr << "readLas: failed to open stream" << std::endl;
    return false;
  }

  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  liblas::Header const &header = reader.GetHeader();

  Eigen::Vector3d offset(header.GetOffsetX(), header.GetOffsetY(), header.GetOffsetZ());
  if (offset_to_remove)
  {
    *offset_to_remove = offset;
    std::cout << "offset to remove: " << offset.transpose() << std::endl;
  }

  const size_t number_of_points = header.GetPointRecordsCount();

  bool using_time = (header.GetDataFormatId() & 1) > 0;
  bool using_colour = (header.GetDataFormatId() & 2) > 0;
  if (!using_time)
  {
    std::cerr << "No timestamps found on laz file, these are required" << std::endl;
    return false;
  }

  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);
  const size_t num_chunks = (number_of_points + (chunk_size - 1)) / chunk_size;
  chunk_size = std::min(number_of_points, chunk_size);
  progress.begin("read and process", num_chunks);

  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends;
  std::vector<double> times;
  std::vector<RGBA> colours;
  std::vector<uint8_t> intensities;
  starts.reserve(chunk_size);
  ends.reserve(chunk_size);
  times.reserve(chunk_size);
  intensities.reserve(chunk_size);
  colours.reserve(chunk_size);

  num_bounded = 0;
  for (unsigned int i = 0; i < number_of_points; i++)
  {
    reader.ReadNextPoint();
    liblas::Point point = reader.GetPoint();

    Eigen::Vector3d position;
    position[0] = point.GetX();
    position[1] = point.GetY();
    position[2] = point.GetZ();

    ends.push_back(position);
    starts.push_back(position);  // equal to position for laz files, as we do not store the start points

    if (using_colour)
    {
      liblas::Color colour = point.GetColor();
      RGBA col;
      col.red = static_cast<uint8_t>(colour.GetRed());
      col.green = static_cast<uint8_t>(colour.GetGreen());
      col.blue = static_cast<uint8_t>(colour.GetBlue());
      colours.push_back(col);
    }
    times.push_back(point.GetTime());

    const double point_int = point.GetIntensity();
    const double normalised_intensity = (255.0 * point_int) / max_intensity;
    const uint8_t intensity = static_cast<uint8_t>(std::min(normalised_intensity, 255.0));
    if (intensity > 0)
      num_bounded++;
    intensities.push_back(intensity);

    if (ends.size() == chunk_size || i == number_of_points - 1)
    {
      if (colours.size() == 0)
      {
        colourByTime(times, colours);
      }
      for (int i = 0; i < (int)colours.size(); i++)  // add intensity into alhpa channel
        colours[i].alpha = intensities[i];
      apply(starts, ends, times, colours);
      starts.clear();
      ends.clear();
      times.clear();
      colours.clear();
      intensities.clear();
      progress.increment();
    }
  }

  progress.end();
  progress_thread.requestQuit();
  progress_thread.join();

  std::cout << "loaded " << file_name << " with " << number_of_points << " points" << std::endl;
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(apply);
  RAYLIB_UNUSED(num_bounded);
  RAYLIB_UNUSED(chunk_size);
  RAYLIB_UNUSED(max_intensity);
  std::cerr << "readLas: cannot read file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
             std::vector<RGBA> &colours, double max_intensity, Eigen::Vector3d *offset_to_remove)
{
  std::vector<Eigen::Vector3d> starts;  // dummy as lax just reads in point clouds, not ray clouds
  auto apply = [&](std::vector<Eigen::Vector3d> &start_points, std::vector<Eigen::Vector3d> &end_points,
                   std::vector<double> &time_points, std::vector<RGBA> &colour_values) 
  {
    starts.insert(starts.end(), start_points.begin(), start_points.end());
    positions.insert(positions.end(), end_points.begin(), end_points.end());
    times.insert(times.end(), time_points.begin(), time_points.end());
    colours.insert(colours.end(), colour_values.begin(), colour_values.end());    
  };
  size_t num_bounded;
  bool success =
    readLas(file_name, apply, num_bounded, max_intensity, offset_to_remove, std::numeric_limits<size_t>::max());
  if (num_bounded == 0)
  {
    std::cout << "warning: all laz file intensities are 0, which would make all rays unbounded. Setting them to 1."
              << std::endl;
    for (auto &c : colours) c.alpha = 255;
  }
  return success;
}

bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points,
                            const std::vector<double> &times, const std::vector<RGBA> &colours)
{
#if RAYLIB_WITH_LAS
  std::cout << "saving LAZ file" << std::endl;

  liblas::Header header;
  header.SetDataFormatId(liblas::ePointFormat1);  // Time only

  if (file_name.find(".laz") != std::string::npos)
    header.SetCompressed(true);

  std::cout << "Saving points to " << file_name << std::endl;

  std::ofstream ofs;
  ofs.open(file_name.c_str(), std::ios::out | std::ios::binary);
  if (ofs.fail())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing." << std::endl;
    return false;
  }

  const double scale = 1e-4;
  header.SetScale(scale, scale, scale);

  liblas::Writer writer(ofs, header);
  liblas::Point point(&header);
  point.SetHeader(&header);  // TODO HACK Version 1.7.0 does not correctly resize the data. Commit
                             // 6e8657336ba445fcec3c9e70c2ebcd2e25af40b9 (1.8.0 3 July fixes it)
  for (unsigned int i = 0; i < points.size(); i++)
  {
    point.SetCoordinates(points[i][0], points[i][1], points[i][2]);
    point.SetIntensity(colours[i].alpha);
    if (!times.empty())
      point.SetTime(times[i]);
    writer.WritePoint(point);
  }
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

#if RAYLIB_WITH_LAS
LasWriter::LasWriter(const std::string &file_name)
  : file_name_(file_name)
{
  header_.SetDataFormatId(liblas::ePointFormat1);  // Time only
  if (file_name_.find(".laz") != std::string::npos)
    header_.SetCompressed(true);

  std::cout << "Saving points to " << file_name_ << std::endl;
  out_.open(file_name_.c_str(), std::ios::out | std::ios::binary);
  if (out_.fail())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing." << std::endl;
    return;
  }
  const double scale = 1e-4;
  header_.SetScale(scale, scale, scale);
  writer_ = new liblas::Writer(out_, header_);
}
#else   // RAYLIB_WITH_LAS
LasWriter::LasWriter(const std::string &file_name)
  : file_name_(file_name)
{
  RAYLIB_UNUSED(file_name);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
}
#endif  // RAYLIB_WITH_LAS

LasWriter::~LasWriter()
{
#if RAYLIB_WITH_LAS
  delete writer_;
#else
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
#endif
}

bool LasWriter::writeChunk(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                           const std::vector<RGBA> &colours)
{
#if RAYLIB_WITH_LAS
  if (points.size() == 0)
  {
    return true;  // this is acceptable behaviour. It avoids calling function checking for emptiness each time
  }
  if (out_.fail())
  {
    std::cerr << "Error: cannot open " << file_name_ << " for writing." << std::endl;
    return false;
  }
  liblas::Point point(&header_);
  point.SetHeader(&header_);  // TODO HACK Version 1.7.0 does not correctly resize the data. Commit
                              // 6e8657336ba445fcec3c9e70c2ebcd2e25af40b9 (1.8.0 3 July fixes it)
  for (unsigned int i = 0; i < points.size(); i++)
  {
    point.SetCoordinates(points[i][0], points[i][1], points[i][2]);
    point.SetIntensity(colours[i].alpha);
    if (!times.empty())
      point.SetTime(times[i]);
    writer_->WritePoint(point);
  }
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

}  // namespace ray
