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
#include <laszip/laszip_api.h>
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

  laszip_POINTER reader;
  if (laszip_create(&reader))
  {
    std::cerr << "readLas: failed to create LASzip reader" << std::endl;
    return false;
  }

  laszip_BOOL is_compressed;
  if (laszip_open_reader(reader, file_name.c_str(), &is_compressed))
  {
    laszip_CHAR *error;
    laszip_get_error(reader, &error);
    std::cerr << "readLas: failed to open stream: " << error << std::endl;
    laszip_destroy(reader);
    return false;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(reader, &header);

  Eigen::Vector3d offset(header->x_offset, header->y_offset, header->z_offset);
  if (offset_to_remove)
  {
    *offset_to_remove = offset;
    std::cout << "offset to remove: " << offset.transpose() << std::endl;
  }

  // LAS 1.4 uses a 64-bit point count; legacy uses the 32-bit field
  const size_t number_of_points =
    (header->version_minor >= 4 && header->extended_number_of_point_records > 0)
      ? static_cast<size_t>(header->extended_number_of_point_records)
      : static_cast<size_t>(header->number_of_point_records);

  const uint8_t format = header->point_data_format;
  // Formats 1,3,4,5 have GPS time in LAS 1.0-1.3; formats 6-10 always have GPS time (LAS 1.4)
  const bool using_time = (format == 1 || format == 3 || format == 4 || format == 5 || format >= 6);
  // Formats 2,3,5 have RGB in LAS 1.0-1.3; formats 7,8,10 have RGB in LAS 1.4
  const bool using_colour = (format == 2 || format == 3 || format == 5 || format == 7 || format == 8 || format == 10);

  if (!using_time)
  {
    std::cerr << "No timestamps found on laz file, these are required" << std::endl;
    laszip_close_reader(reader);
    laszip_destroy(reader);
    return false;
  }

  laszip_point_struct *point;
  laszip_get_point_pointer(reader, &point);

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
  for (size_t i = 0; i < number_of_points; i++)
  {
    if (laszip_read_point(reader))
    {
      laszip_CHAR *error;
      laszip_get_error(reader, &error);
      std::cerr << "readLas: error reading point " << i << ": " << error << std::endl;
      break;
    }

    laszip_F64 coords[3];
    laszip_get_coordinates(reader, coords);
    Eigen::Vector3d position(coords[0], coords[1], coords[2]);

    ends.push_back(position);
    starts.push_back(position);  // equal to position for laz files, as we do not store the start points

    if (using_colour)
    {
      RGBA col;
      col.red = static_cast<uint8_t>(point->rgb[0]);
      col.green = static_cast<uint8_t>(point->rgb[1]);
      col.blue = static_cast<uint8_t>(point->rgb[2]);
      colours.push_back(col);
    }
    times.push_back(point->gps_time);

    const double point_int = point->intensity;
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
      for (size_t j = 0; j < colours.size(); j++)  // add intensity into alpha channel
        colours[j].alpha = intensities[j];
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

  laszip_close_reader(reader);
  laszip_destroy(reader);

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

  laszip_POINTER writer;
  if (laszip_create(&writer))
  {
    std::cerr << "writeLas: failed to create LASzip writer" << std::endl;
    return false;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(writer, &header);

  header->version_major = 1;
  header->version_minor = 2;
  header->point_data_format = 1;  // GPS time only
  const double scale = 1e-4;
  header->x_scale_factor = scale;
  header->y_scale_factor = scale;
  header->z_scale_factor = scale;
  header->x_offset = 0.0;
  header->y_offset = 0.0;
  header->z_offset = 0.0;
  header->number_of_point_records = static_cast<laszip_U32>(points.size());

  const bool is_laz = file_name.find(".laz") != std::string::npos;
  std::cout << "Saving points to " << file_name << std::endl;

  if (laszip_open_writer(writer, file_name.c_str(), is_laz ? 1 : 0))
  {
    laszip_CHAR *error;
    laszip_get_error(writer, &error);
    std::cerr << "writeLas: failed to open file for writing: " << error << std::endl;
    laszip_destroy(writer);
    return false;
  }

  laszip_point_struct *point;
  laszip_get_point_pointer(writer, &point);

  for (size_t i = 0; i < points.size(); i++)
  {
    laszip_F64 coords[3] = { points[i][0], points[i][1], points[i][2] };
    laszip_set_coordinates(writer, coords);
    point->intensity = colours[i].alpha;
    if (!times.empty())
      point->gps_time = times[i];
    laszip_write_point(writer);
  }

  laszip_update_inventory(writer);
  laszip_close_writer(writer);
  laszip_destroy(writer);
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
  , writer_handle_(nullptr)
  , header_(nullptr)
  , point_(nullptr)
{
  if (laszip_create(&writer_handle_))
  {
    std::cerr << "LasWriter: failed to create LASzip writer" << std::endl;
    writer_handle_ = nullptr;
    return;
  }

  laszip_get_header_pointer(writer_handle_, &header_);

  header_->version_major = 1;
  header_->version_minor = 2;
  header_->point_data_format = 1;  // GPS time only
  const double scale = 1e-4;
  header_->x_scale_factor = scale;
  header_->y_scale_factor = scale;
  header_->z_scale_factor = scale;
  header_->x_offset = 0.0;
  header_->y_offset = 0.0;
  header_->z_offset = 0.0;

  const bool is_laz = file_name_.find(".laz") != std::string::npos;
  std::cout << "Saving points to " << file_name_ << std::endl;

  if (laszip_open_writer(writer_handle_, file_name_.c_str(), is_laz ? 1 : 0))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasWriter: failed to open file for writing: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  laszip_get_point_pointer(writer_handle_, &point_);
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
  if (writer_handle_)
  {
    laszip_update_inventory(writer_handle_);
    laszip_close_writer(writer_handle_);
    laszip_destroy(writer_handle_);
  }
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
  if (!writer_handle_ || !point_)
  {
    std::cerr << "Error: cannot open " << file_name_ << " for writing." << std::endl;
    return false;
  }
  for (size_t i = 0; i < points.size(); i++)
  {
    laszip_F64 coords[3] = { points[i][0], points[i][1], points[i][2] };
    laszip_set_coordinates(writer_handle_, coords);
    point_->intensity = colours[i].alpha;
    if (!times.empty())
      point_->gps_time = times[i];
    laszip_write_point(writer_handle_);
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
