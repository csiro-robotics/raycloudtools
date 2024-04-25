// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230

// Author: Tim Devereux
#include "rayriegl.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raypose.h"
#include "rayunused.h"


#if RAYLIB_WITH_RIEGL
#include <cmath>
#include <iostream>
// #include <riegl/rdb.hpp>
// #include <riegl/rdb/default.hpp>
#include <riegl/scanlib.hpp>
#include <tuple>
#endif  // RAYLIB_WITH_RIEGL


struct RieglPointData
{
  std::vector<float> start_x, start_y, start_z;
  std::vector<float> end_x, end_y, end_z;
  std::vector<float> reflectance;
  std::vector<float> time;
};

std::tuple<float, float, float> normalizeVector(float Vx, float Vy, float Vz)
{
  float magnitude = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
  return std::make_tuple(Vx / magnitude, Vy / magnitude, Vz / magnitude);
}

std::tuple<float, float, float> getPositionAtDistance(float X0, float Y0, float Z0, float Vx, float Vy, float Vz,
                                                      float distance)
{
  double Vx_normalized, Vy_normalized, Vz_normalized;
  std::tie(Vx_normalized, Vy_normalized, Vz_normalized) = normalizeVector(Vx, Vy, Vz);
  float X = X0 + distance * Vx_normalized;
  float Y = Y0 + distance * Vy_normalized;
  float Z = Z0 + distance * Vz_normalized;
  return std::make_tuple(X, Y, Z);
}

class RieglReader : public scanlib::pointcloud
{
  RieglPointData &rxp_data;

public:
  RieglReader(RieglPointData &rxp_data)
    : scanlib::pointcloud(false)
    , rxp_data(rxp_data)
  {}

protected:
  void on_shot_end()
  {
    if (target_count == 0 && beam_direction[2] > 0 && beam_direction[2] < 0.866) //remove points below the scanner and remove the "buffer" points at the top of the scan lines
    {
      float X, Y, Z;
      std::tie(X, Y, Z) = getPositionAtDistance(beam_origin[0], beam_origin[1], beam_origin[2], beam_direction[0],
                                                beam_direction[1], beam_direction[2], 1000);
      rxp_data.end_x.push_back(X);
      rxp_data.end_y.push_back(Y);
      rxp_data.end_z.push_back(Z);
      rxp_data.start_x.push_back(beam_origin[0]);
      rxp_data.start_y.push_back(beam_origin[1]);
      rxp_data.start_z.push_back(beam_origin[2]);
      rxp_data.reflectance.push_back(0);
      rxp_data.time.push_back(time);
    }
  }
};

namespace ray
{
bool readRXP(const std::string &file_name,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply,
             size_t &num_bounded, double max_intensity, std::vector<double> pose_transformation, size_t chunk_size)
{
  std::cout << "readRiegl: " << file_name << std::endl;

  std::ifstream ifs;
  ifs.open(file_name.c_str(), std::ios::in | std::ios::binary);

  if (ifs.fail())
  {
    std::cerr << "readRiegl: failed to open stream" << std::endl;
    return false;
  }

  RieglPointData riegl_data;

  std::shared_ptr<scanlib::basic_rconnection> rc;
  rc = scanlib::basic_rconnection::create(file_name);
  rc->open();
  scanlib::decoder_rxpmarker dec(rc);
  RieglReader imp(riegl_data);
  scanlib::buffer buf;
  rc->close();
  for (dec.get(buf); !dec.eoi(); dec.get(buf))
  {
    imp.dispatch(buf.begin(), buf.end());
  }

  std::size_t number_of_points = riegl_data.start_x.size();

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
    Eigen::Vector3d start;
    Eigen::Vector3d end;

    start[0] = ((riegl_data.end_x[i] * pose_transformation[0]) + (riegl_data.end_y[i] * pose_transformation[1]) +
                (riegl_data.end_z[i] * pose_transformation[2])) +
               pose_transformation[3];
    start[1] = ((riegl_data.end_x[i] * pose_transformation[4]) + (riegl_data.end_y[i] * pose_transformation[5]) +
                (riegl_data.end_z[i] * pose_transformation[6])) +
               pose_transformation[7];
    start[2] = ((riegl_data.end_x[i] * pose_transformation[8]) + (riegl_data.end_y[i] * pose_transformation[9]) +
                (riegl_data.end_z[i] * pose_transformation[10])) +
               pose_transformation[11];

    end[0] = ((riegl_data.end_x[i] * pose_transformation[0]) + (riegl_data.end_y[i] * pose_transformation[1]) +
              (riegl_data.end_z[i] * pose_transformation[2])) +
             pose_transformation[3];
    end[1] = ((riegl_data.end_x[i] * pose_transformation[4]) + (riegl_data.end_y[i] * pose_transformation[5]) +
              (riegl_data.end_z[i] * pose_transformation[6])) +
             pose_transformation[7];
    end[2] = ((riegl_data.end_x[i] * pose_transformation[8]) + (riegl_data.end_y[i] * pose_transformation[9]) +
              (riegl_data.end_z[i] * pose_transformation[10])) +
             pose_transformation[11];

    ends.push_back(end);
    starts.push_back(start);
    times.push_back(riegl_data.time[i]);

    const double point_int = riegl_data.reflectance[i];
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
}


// bool readRDBX(const std::string &file_name,
//               std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
//                                  std::vector<double> &times, std::vector<RGBA> &colours)>
//                 apply,
//               size_t &num_bounded, double max_intensity, size_t chunk_size)
// {
  
//   // New RDB library context
//   riegl::rdb::Context context;

//   // Access existing database
//   riegl::rdb::Pointcloud rdb(context);
//   riegl::rdb::pointcloud::OpenSettings settings;
//   rdb.open(file_name, settings);

//   // Prepare point attribute buffers
//   static const uint32_t BUFFER_SIZE = 10000;
//   std::vector<uint64_t> bufferPointNumber(BUFFER_SIZE);
//   std::vector<std::array<double, 3>> bufferCoordinates(BUFFER_SIZE);
//   std::vector<float> bufferReflectance(BUFFER_SIZE);
//   std::vector<std::array<uint8_t, 4>> bufferTimestamp(BUFFER_SIZE);
//   std::vector<std::array<uint8_t, 4>> bufferTrueColor(BUFFER_SIZE);

//   // Start new select query to read...
//   riegl::rdb::pointcloud::QuerySelect select = rdb.select();

//   // Tell select query where to store the data
//   using namespace riegl::rdb::pointcloud;
//   select.bindBuffer(RDB_RIEGL_ID, bufferPointNumber);
//   select.bindBuffer(RDB_RIEGL_XYZ, bufferCoordinates);
//   select.bindBuffer(RDB_RIEGL_REFLECTANCE, bufferReflectance);
//   // select.bindBuffer(RDB_RIEGL_TIMESTAMP, bufferTimestamp);
//   // select.bindBuffer(RDB_RIEGL_RGBA, bufferTrueColor);

//   // ray::Progress progress;
//   // ray::ProgressThread progress_thread(progress);
//   // progress.begin("read and process", chunk_size);

//   std::vector<Eigen::Vector3d> starts;
//   std::vector<Eigen::Vector3d> ends;
//   std::vector<double> times;
//   std::vector<RGBA> colours;
//   std::vector<uint8_t> intensities;
//   starts.reserve(chunk_size);
//   ends.reserve(chunk_size);
//   times.reserve(chunk_size);
//   // intensities.reserve(chunk_size);
//   // colours.reserve(chunk_size);

//   num_bounded = 0;
//   unsigned int number_of_points = bufferPointNumber[-1];

//   // Read and process all points block-wise
//   while (const uint32_t count = select.next(BUFFER_SIZE))
//   {
//     // Print points to output stream
//     for (uint32_t i = 0; i < count; i++)
//     {
//       Eigen::Vector3d end;
//       end[0] = bufferCoordinates[i][0];
//       end[1] = bufferCoordinates[i][1];
//       end[2] = bufferCoordinates[i][2];

//       // std::vector<RGBA> colour;
//       // colour << int(bufferTrueColor[i][0]), int(bufferTrueColor[i][0]), int(bufferTrueColor[i][0]);
//       // ends.push_back(end);
//       // times.push_back(bufferTimestamp[i]);
//       // colours.push_back(colour);

//       const double point_int = bufferReflectance[i];
//       const double normalised_intensity = (255.0 * point_int) / max_intensity;
//       const uint8_t intensity = static_cast<uint8_t>(std::min(normalised_intensity, 255.0));
//       if (point_int > 0)
//         num_bounded++;
//       intensities.push_back(point_int);


//       if (ends.size() == BUFFER_SIZE || i == number_of_points - 1)
//       {
//         if (colours.size() == 0)
//         {
//           colourByTime(times, colours);
//         }
//         for (int i = 0; i < (int)colours.size(); i++)  // add intensity into alhpa channel
//           colours[i].alpha = intensities[i];
//         apply(starts, ends, times, colours);
//         starts.clear();
//         ends.clear();
//         times.clear();
//         colours.clear();
//         intensities.clear();
//         // progress.increment();
//       }
//     }
//   }
//   // progress.end();
//   // progress_thread.requestQuit();
//   // progress_thread.join();

//   std::cout << "loaded " << file_name << " with " << number_of_points << " points" << std::endl;

//   return true;
// }


}  // namespace ray