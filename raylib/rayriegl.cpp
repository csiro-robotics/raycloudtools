// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
// Author: Tim Devereux

#include "rayriegl.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raypose.h"
#include "rayunused.h"

namespace ray
{

#if RAYLIB_WITH_RIEGL
#include <cmath>
#include <iostream>
#include <memory>
#include <riegl/scanlib.hpp>
#include <stdexcept>
#include <tuple>

namespace
{

struct RieglPointData
{
  std::vector<float> start_x, start_y, start_z;
  std::vector<float> end_x, end_y, end_z;
  std::vector<float> time;

  size_t size() const { return start_x.size(); }

  void reserve(size_t capacity)
  {
    start_x.reserve(capacity);
    start_y.reserve(capacity);
    start_z.reserve(capacity);
    end_x.reserve(capacity);
    end_y.reserve(capacity);
    end_z.reserve(capacity);
    time.reserve(capacity);
  }

  void clear()
  {
    start_x.clear();
    start_y.clear();
    start_z.clear();
    end_x.clear();
    end_y.clear();
    end_z.clear();
    time.clear();
  }
};

class Vector3D
{
public:
  static std::tuple<float, float, float> normalize(float x, float y, float z)
  {
    float magnitude = std::sqrt(x * x + y * y + z * z);
    return std::make_tuple(x / magnitude, y / magnitude, z / magnitude);
  }

  static std::tuple<float, float, float> getPositionAtDistance(float x0, float y0, float z0, float vx, float vy,
                                                               float vz, float distance)
  {
    auto [vx_norm, vy_norm, vz_norm] = normalize(vx, vy, vz);
    return std::make_tuple(x0 + distance * vx_norm, y0 + distance * vy_norm, z0 + distance * vz_norm);
  }
};
class RieglReader : public scanlib::pointcloud {
public:
    explicit RieglReader(RieglPointData& rxp_data)
        : scanlib::pointcloud(false)
        , rxp_data_(rxp_data) 
    {}

protected:
    void on_shot_end() override {
        static constexpr float MIN_Z_DIRECTION = 0.0f;
        static constexpr float MAX_Z_DIRECTION = 0.866f;
        
        // Skip points below scanner and buffer points at top of scan lines
        if (target_count != 0 || 
            // beam_direction[2] <= MIN_Z_DIRECTION ||
            beam_direction[2] >= MAX_Z_DIRECTION) {
            return;
        }

        auto [end_x, end_y, end_z] = Vector3D::getPositionAtDistance(
            beam_origin[0], beam_origin[1], beam_origin[2],
            beam_direction[0], beam_direction[1], beam_direction[2],
            1000.0f
        );

        rxp_data_.end_x.push_back(end_x);
        rxp_data_.end_y.push_back(end_y);
        rxp_data_.end_z.push_back(end_z);
        rxp_data_.start_x.push_back(beam_origin[0]);
        rxp_data_.start_y.push_back(beam_origin[1]);
        rxp_data_.start_z.push_back(beam_origin[2]);
        rxp_data_.time.push_back(time);
    }

private:
    RieglPointData& rxp_data_;
};

class RieglDataProcessor
{
public:
  static Eigen::Vector3d transformPoint(const std::vector<double> &pose_transform, float x, float y, float z)
  {
    if (pose_transform.empty())
    {
      return Eigen::Vector3d(x, y, z);
    }

    Eigen::Vector3d result;
    for (int i = 0; i < 3; ++i)
    {
      result[i] = (x * pose_transform[i * 4] + y * pose_transform[i * 4 + 1] + z * pose_transform[i * 4 + 2]) +
                  pose_transform[i * 4 + 3];
    }
    return result;
  }
};

}  // anonymous namespace

bool readRXP(const std::string& file_name,
             std::function<void(std::vector<Eigen::Vector3d>& starts,
                              std::vector<Eigen::Vector3d>& ends,
                              std::vector<double>& times,
                              std::vector<RGBA>& colours)> apply,
             size_t& num_bounded,
             [[maybe_unused]] double max_intensity,
             std::vector<double> pose_transformation,
             size_t chunk_size)
{
    try {
        std::cout << "Reading Riegl file: " << file_name << std::endl;
        
        auto rc = scanlib::basic_rconnection::create(file_name);
        if (!rc) {
            throw std::runtime_error("Failed to create connection");
        }
        
        RieglPointData riegl_data;
        rc->open();
        scanlib::decoder_rxpmarker dec(rc);
        RieglReader reader(riegl_data);
        scanlib::buffer buf;
        
        for (dec.get(buf); !dec.eoi(); dec.get(buf)) {
            reader.dispatch(buf.begin(), buf.end());
        }
        rc->close();

        const size_t number_of_points = riegl_data.size();
        if (number_of_points == 0) {
            throw std::runtime_error("No points read from file");
        }

        ray::Progress progress;
        ray::ProgressThread progress_thread(progress);
        const size_t num_chunks = (number_of_points + (chunk_size - 1)) / chunk_size;
        chunk_size = std::min(number_of_points, chunk_size);
        progress.begin("Processing points", num_chunks);

        std::vector<Eigen::Vector3d> starts;
        std::vector<Eigen::Vector3d> ends;
        std::vector<double> times;
        std::vector<RGBA> colours;

        starts.reserve(chunk_size);
        ends.reserve(chunk_size);
        times.reserve(chunk_size);
        colours.reserve(chunk_size);

        num_bounded = 0;  // Set to 0 since we're not tracking intensity bounds

        for (size_t i = 0; i < number_of_points; ++i) {
            auto start = RieglDataProcessor::transformPoint(
                pose_transformation,
                riegl_data.start_x[i],
                riegl_data.start_y[i],
                riegl_data.start_z[i]
            );

            auto end = RieglDataProcessor::transformPoint(
                pose_transformation,
                riegl_data.end_x[i],
                riegl_data.end_y[i],
                riegl_data.end_z[i]
            );

            ends.push_back(end);
            starts.push_back(start);
            times.push_back(riegl_data.time[i]);

            if (ends.size() == chunk_size || i == number_of_points - 1) {
                if (colours.empty()) {
                    colourByTime(times, colours);
                    // Set all intensities to 0
                    colours.resize(times.size());
                    for (auto& colour : colours) {
                        colour.alpha = u_int8_t(0);
                    }
                }
                
                apply(starts, ends, times, colours);
                
                starts.clear();
                ends.clear();
                times.clear();
                colours.clear();
                
                progress.increment();
            }
        }

        progress.end();
        progress_thread.requestQuit();
        progress_thread.join();

        std::cout << "Successfully loaded " << number_of_points << " points from " 
                  << file_name << std::endl;
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error reading RXP file: " << e.what() << std::endl;
        return false;
    }
}
#else  // !RAYLIB_WITH_RIEGL

bool readRXP(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
             std::vector<RGBA> &colours, double max_intensity, std::vector<double> pose_transformation)
{
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(positions);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(pose_transformation);
  std::cerr << "RIEGL support is not enabled in this build" << std::endl;
  return false;
}

bool readRXP(const std::string &file_name,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply,
             size_t &num_bounded, double max_intensity, std::vector<double> pose_transformation, size_t chunk_size)
{
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(apply);
  RAYLIB_UNUSED(num_bounded);
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(pose_transformation);
  RAYLIB_UNUSED(chunk_size);
  std::cerr << "RIEGL support is not enabled in this build" << std::endl;
  return false;
}

#endif  // RAYLIB_WITH_RIEGL

}  // namespace ray