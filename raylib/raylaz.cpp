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
#include "laszip_api.h"
#endif  // RAYLIB_WITH_LAS

namespace ray
{
/// Write header fields to std::cout.
/// For debugging.
void printHeader(const std::string &prefix, const laszip_header_struct *const header)
{
  std::cout << prefix << "file_source_ID: " << header->file_source_ID << "\n";
  std::cout << prefix << "global_encoding: " << header->global_encoding << "\n";
  std::cout << prefix << "project_ID_GUID_data_1: " << header->project_ID_GUID_data_1 << "\n";
  std::cout << prefix << "project_ID_GUID_data_2: " << header->project_ID_GUID_data_2 << "\n";
  std::cout << prefix << "project_ID_GUID_data_3: " << header->project_ID_GUID_data_3 << "\n";
  std::cout << prefix << "project_ID_GUID_data_4: " << header->project_ID_GUID_data_4 << "\n";
  std::cout << prefix << "version_major: " << static_cast<int>(header->version_major) << "\n";
  std::cout << prefix << "version_minor: " << static_cast<int>(header->version_minor) << "\n";
  std::cout << prefix << "system_identifier: " << header->system_identifier << "\n";
  std::cout << prefix << "generating_software: " << header->generating_software << "\n";
  std::cout << prefix << "file_creation_day: " << header->file_creation_day << "\n";
  std::cout << prefix << "file_creation_year: " << header->file_creation_year << "\n";
  std::cout << prefix << "header_size: " << header->header_size << "\n";
  std::cout << prefix << "offset_to_point_data: " << header->offset_to_point_data << "\n";
  std::cout << prefix << "number_of_variable_length_records: " << header->number_of_variable_length_records << "\n";
  std::cout << prefix << "point_data_format: " << static_cast<int>(header->point_data_format) << "\n";
  std::cout << prefix << "point_data_record_length: " << header->point_data_record_length << "\n";
  std::cout << prefix << "number_of_point_records: " << header->number_of_point_records << "\n";
  for (size_t index = 0; index < 5; ++index)
  {
    std::cout << prefix << "number_of_points_by_return[" << index << "]: " << header->number_of_points_by_return[index]
              << "\n";
  }
  std::cout << prefix << "x_scale_factor: " << header->x_scale_factor << "\n";
  std::cout << prefix << "y_scale_factor: " << header->y_scale_factor << "\n";
  std::cout << prefix << "z_scale_factor: " << header->z_scale_factor << "\n";
  std::cout << prefix << "x_offset: " << header->x_offset << "\n";
  std::cout << prefix << "y_offset: " << header->y_offset << "\n";
  std::cout << prefix << "z_offset: " << header->z_offset << "\n";
  std::cout << prefix << "max_x: " << header->max_x << "\n";
  std::cout << prefix << "min_x: " << header->min_x << "\n";
  std::cout << prefix << "max_y: " << header->max_y << "\n";
  std::cout << prefix << "min_y: " << header->min_y << "\n";
  std::cout << prefix << "max_z: " << header->max_z << "\n";
  std::cout << prefix << "min_z: " << header->min_z << "\n";
  std::cout << prefix << "start_of_waveform_data_packet_record: " << header->start_of_waveform_data_packet_record
            << "\n";
  std::cout << prefix << "start_of_first_extended_variable_length_record: "
            << header->start_of_first_extended_variable_length_record << "\n";
  std::cout << prefix
            << "number_of_extended_variable_length_records: " << header->number_of_extended_variable_length_records
            << "\n";
  std::cout << prefix << "extended_number_of_point_records: " << header->extended_number_of_point_records << "\n";
  for (size_t index = 0; index < 15; ++index)
  {
    std::cout << prefix << "extended_number_of_points_by_return[" << index
              << "]: " << header->extended_number_of_points_by_return[index] << "\n";
  }
  // std::cout << prefix << "max_gps_time: " << header->max_gps_time << "\n";  // Requires LASzip >= 3.5.0.
  // std::cout << prefix << "min_gps_time: " << header->min_gps_time << "\n";  // Requires LASzip >= 3.5.0.
  // std::cout << prefix << "time_offset: " << header->time_offset << "\n";    // Requires LASzip >= 3.5.0.
  std::cout << prefix << "user_data_in_header_size: " << header->user_data_in_header_size << "\n";
  std::cout << prefix << "user_data_after_header_size: " << header->user_data_after_header_size << "\n";
}

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

  // printHeader("read_header.", header);

  Eigen::Vector3d offset(header->x_offset, header->y_offset, header->z_offset);
  if (offset_to_remove)
  {
    *offset_to_remove = offset;
    std::cout << "offset to remove: " << offset.transpose() << std::endl;
  }

  // LAS 1.4 uses a 64-bit point count; legacy uses the 32-bit field
  const size_t number_of_points = (header->version_minor >= 4 && header->extended_number_of_point_records > 0) ?
                                    static_cast<size_t>(header->extended_number_of_point_records) :
                                    static_cast<size_t>(header->number_of_point_records);

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
      col.red = static_cast<uint8_t>(point->rgb[0] & 0x00FF);
      col.green = static_cast<uint8_t>(point->rgb[1] & 0x00FF);
      col.blue = static_cast<uint8_t>(point->rgb[2] & 0x00FF);
      col.alpha = static_cast<uint8_t>(point->rgb[3] & 0x00FF);
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
  RAYLIB_UNUSED(offset_to_remove);
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(apply);
  RAYLIB_UNUSED(num_bounded);
  RAYLIB_UNUSED(chunk_size);
  RAYLIB_UNUSED(max_intensity);
  std::cerr << "readLas: cannot read file as RAYLIB_WITH_LAS not enabled. "
            << "Enable using: cmake .. -DRAYLIB_WITH_LAS:BOOL=ON" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
             std::vector<RGBA> &colours, double max_intensity, Eigen::Vector3d *offset_to_remove)
{
  std::vector<Eigen::Vector3d> starts;  // dummy as lax just reads in point clouds, not ray clouds
  auto apply = [&](std::vector<Eigen::Vector3d> &start_points, std::vector<Eigen::Vector3d> &end_points,
                   std::vector<double> &time_points, std::vector<RGBA> &colour_values) {
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
  auto writer = LasWriter(file_name);
  return writer.writeChunk(points, times, colours);
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as RAYLIB_WITH_LAS not enabled. "
            << "Enable using: cmake .. -DRAYLIB_WITH_LAS:BOOL=ON" << std::endl;
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

  // https://www.lvermgeo.sachsen-anhalt.de/datei/anzeigen/id/624579,501/asprs_las_format_v12.pdf
  // https://downloads.rapidlasso.de/doc/LAZ_Specification_1.4_R1.pdf

  laszip_set_point_type_and_size(writer_handle_, 3, 34);  // Point10 + GPSTime11 + RGB12

  const double scale = 1e-4;
  header_->x_scale_factor = scale;
  header_->y_scale_factor = scale;
  header_->z_scale_factor = scale;
  header_->x_offset = 0.0;
  header_->y_offset = 0.0;
  header_->z_offset = 0.0;

  // printHeader("write_header.", header_);

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
    if (laszip_close_writer(writer_handle_))
    {
      std::cerr << "Error: laszip_close_writer failed.\n";
    }
    if (laszip_destroy(writer_handle_))
    {
      std::cerr << "Error: laszip_destroy failed.\n";
    }
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
    // This is acceptable behaviour. It avoids calling function checking for
    // emptiness each time.
    return true;
  }
  if (!writer_handle_ || !point_)
  {
    std::cerr << "Error: cannot open " << file_name_ << " for writing." << std::endl;
    return false;
  }
  for (size_t i = 0; i < points.size(); i++)
  {
    const laszip_F64 coords[3] = { points[i][0], points[i][1], points[i][2] };
    laszip_set_coordinates(writer_handle_, coords);
    point_->intensity = colours[i].alpha;
    if (!times.empty())
    {
      point_->gps_time = times[i];
    }
    if (!colours.empty())
    {
      point_->rgb[0] = static_cast<laszip_U16>(colours[i].red);
      point_->rgb[1] = static_cast<laszip_U16>(colours[i].green);
      point_->rgb[2] = static_cast<laszip_U16>(colours[i].blue);
      point_->rgb[3] = static_cast<laszip_U16>(colours[i].alpha);
    }
    if (laszip_write_point(writer_handle_))
    {
      std::cerr << "Error: laszip_write_point failed.\n";
      return false;
    }
    if (laszip_update_inventory(writer_handle_))
    {
      std::cerr << "Error: laszip_update_inventory failed.\n";
      return false;
    }
  }
  int64_t wrote_points = -1;
  if (laszip_get_point_count(writer_handle_, &wrote_points))
  {
    std::cerr << "Error: laszip_get_point_count failed.\n";
    return false;
  }
  std::cout << "Wrote " << wrote_points << " points.\n";
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
