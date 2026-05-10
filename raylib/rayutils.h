// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYUTILS_H
#define RAYLIB_RAYUTILS_H

#include "raylib/raylibconfig.h"
#include "rayrandom.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#define VISUALISE_TOOL "QT_QPA_PLATFORM=xcb meshlab" // the first term fixed opening problems on some platforms

namespace ray
{
const double kPi = M_PI;
// while this is an absolute value, it has little effect on results unless point spacing is signficantly less
// than this small value in metres. However, the computation time is better for having a value greater than 0..
const double kNearestNeighbourEpsilon = 0.001;
#define ASSERT(X) assert(X);

inline int runWithMemoryCheck(std::function<bool(int argc, char *argv[])> main_function, int argc, char *argv[])
{
  try
  {
    int result = main_function(argc, argv);
    return result;
  }
  catch (std::bad_alloc const &)  // catch any memory allocation problems in generating large images
  {
    std::cerr << "Error: Not enough memory to process the input file," << std::endl;
    std::cerr << "consider using raydecimate or raysplit grid to operate on a smaller file." << std::endl;
    return 1;
  }  
  catch (std::length_error const &)  // catch any memory allocation problems in generating large images
  {
    std::cerr << "Error: Not enough memory to process the input file," << std::endl;
    std::cerr << "consider using raydecimate or raysplit grid to operate on a smaller file." << std::endl;
    return 1;
  }    
}

inline std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;
  while (getline(ss, item, delim)) result.push_back(item);
  return result;
}

template <class T>
inline const T maxVector(const T &a, const T &b)
{
  return T(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}
template <class T>
inline const T minVector(const T &a, const T &b)
{
  return T(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
}

template <class T>
T clamped(const T &value, const T &min_value, const T &max_value)
{
  return std::max(min_value, std::min(value, max_value));
}

template <typename T>
T sgn(T val)
{
  return (val > T(0)) ? T(1) : T(-1);
}

inline int roundToInt(double x)
{
  if (x >= 0)
    return int(x + 0.5);
  return -int(0.5 - x);
}

/// Uniform distribution within range
inline double random(double min, double max)
{
  return min + ((max - min) * randUniformDouble());
}

class RAYLIB_EXPORT Vector3iLess
{
public:
  bool operator()(const Eigen::Vector3i &a, const Eigen::Vector3i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    return a[2] < b[2];
  }
};

inline void voxelSubsample(const std::vector<Eigen::Vector3d> &points, double voxel_width,
                           std::vector<int64_t> &indices, std::set<Eigen::Vector3i, Vector3iLess> &vox_set)
{
  for (int64_t i = 0; i < (int64_t)points.size(); i++)
  {
    Eigen::Vector3i voxel(int(std::floor(points[i][0] / voxel_width)), int(std::floor(points[i][1] / voxel_width)),
                          int(std::floor(points[i][2] / voxel_width)));
    if (vox_set.insert(voxel).second)
    {
      indices.push_back(i);
    }
  }
}

inline void voxelSubsample(const std::vector<Eigen::Vector3d> &points, double voxel_width,
                           std::vector<int64_t> &indices)
{
  std::set<Eigen::Vector3i, Vector3iLess> vox_set;
  voxelSubsample(points, voxel_width, indices, vox_set);
}

/// Square a value
template <class T>
inline T sqr(const T &val)
{
  return val * val;
}

template <class T>
inline T mean(const std::vector<T> &list)
{
  T result = list[0];
  for (unsigned int i = 1; i < list.size(); i++) result += list[i];
  result /= double(list.size());
  return result;
}

/** Return median of elements in the list
 * When there are an even number of elements it returns the mean of the two medians
 */
template <class T>
inline T median(std::vector<T> list)
{
  typename std::vector<T>::iterator first = list.begin();
  typename std::vector<T>::iterator last = list.end();
  typename std::vector<T>::iterator middle = first + ((last - first) / 2);
  nth_element(first, middle, last);  // can specify comparator as optional 4th arg
  if (list.size() % 2)               // odd
    return *middle;
  else
  {
    typename std::vector<T>::iterator middle2 = middle - 1;
    nth_element(first, middle2, last);
    return (*middle + *middle2) / static_cast<T>(2);
  }
}

/** Returns p'th percentile value in unordered list. e.g. p=50% gives median value, p=0% gives smallest value
 */
template <class T>
T percentile(std::vector<T> list, double p)
{
  typename std::vector<T>::iterator first = list.begin();
  typename std::vector<T>::iterator last = list.end();
  int closest_index = (int)(p * ((double)list.size()) / 100.0);
  typename std::vector<T>::iterator percentile = first + closest_index;
  nth_element(first, percentile, last);  // can specify comparator as optional 4th arg
  return *percentile;
}

// TODO: this function/macro is possibly dangerous due to its casting and manipulation of pointers. Consider
// alternatives. Retrieve a list of a certain component of a vector, e.g. a list of the w values of a vector of
// quaternions: vector<double> wValues = components(quats, quats[0].w);  // vector<Quat> quats Note: this is suboptimal
// in that the list is passed back by value, which requires a copy.
#define components(_list, _component) component_list(_list, _list[0]._component)
template <class U, class T>
inline std::vector<T> componentList(const std::vector<U> &list, const T &component0)
{
  unsigned long int offset = (unsigned long int)&component0 - (unsigned long int)&list[0];
  std::vector<T> sub_list(list.size());
  for (unsigned int i = 0; i < list.size(); ++i) sub_list[i] = *(T *)((char *)&list[i] + offset);
  return sub_list;
}

struct RGBA
{
  RGBA(){}
  RGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : red(r), green(g), blue(b), alpha(a) {}
  uint8_t red;
  uint8_t green;
  uint8_t blue;
  uint8_t alpha;
  static RGBA white(){ return RGBA(255, 255, 255, 255); }
  static RGBA terrain(){ return RGBA(149,105,72, 255); }
  static RGBA treetrunk(){ return RGBA(192,166,141, 255); }
  static RGBA leaves(){ return RGBA(60,102,44, 255); }
};

/// Converts a value from 0 to 1 into a RGBA structure
inline Eigen::Vector3d redGreenBlueGradient(double val)
{
  const Eigen::Vector3d hue_cycle[6] = { Eigen::Vector3d(1.0, 0.0, 0.5), Eigen::Vector3d(1.0, 0.5, 0.0),
                                         Eigen::Vector3d(0.5, 1.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.5),
                                         Eigen::Vector3d(0.0, 0.5, 1.0), Eigen::Vector3d(0.5, 0.0, 1.0) };
  double v = 0.5 + 4.0 * clamped(val, 0.0, 1.0);
  int id = (int)v;
  double blend = v - (double)id;
  return hue_cycle[id] * (1.0 - blend) + hue_cycle[id + 1] * blend;
}

inline void redGreenBlueGradient(const std::vector<double> &values, std::vector<RGBA> &gradient, double min_value,
                                 double max_value, bool replace_alpha)
{
  gradient.resize(values.size());
  for (size_t i = 0; i < values.size(); i++)
  {
    const Eigen::Vector3d col = redGreenBlueGradient((values[i] - min_value) / (max_value - min_value));
    gradient[i].red = uint8_t(255.0 * col[0]);
    gradient[i].green = uint8_t(255.0 * col[1]);
    gradient[i].blue = uint8_t(255.0 * col[2]);
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline Eigen::Vector3d redGreenBlueSpectrum(double value)
{
  const Eigen::Vector3d rgb_cycle[6] = { Eigen::Vector3d(1.0, 0.5, 0.0), Eigen::Vector3d(0.5, 1.0, 0.0),
                                         Eigen::Vector3d(0.0, 1.0, 0.5), Eigen::Vector3d(0.0, 0.5, 1.0),
                                         Eigen::Vector3d(0.5, 0.0, 1.0), Eigen::Vector3d(1.0, 0.0, 0.1) };

  double v = 6.0 * (value - std::floor(value));
  int id = (int)v;
  double blend = v - (double)id;
  return rgb_cycle[id] * (1.0 - blend) + rgb_cycle[(id + 1) % 6] * blend;
}

inline void redGreenBlueSpectrum(const std::vector<double> &values, std::vector<RGBA> &gradient, double wavelength,
                                 bool replace_alpha)
{
  gradient.resize(values.size());
  for (size_t i = 0; i < values.size(); i++)
  {
    const Eigen::Vector3d col = redGreenBlueSpectrum(values[i] / wavelength);
    gradient[i].red = uint8_t(255.0 * col[0]);
    gradient[i].green = uint8_t(255.0 * col[1]);
    gradient[i].blue = uint8_t(255.0 * col[2]);
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline void colourByTime(const std::vector<double> &values, std::vector<RGBA> &gradient, bool replace_alpha = true)
{
  const double colour_repeat_period = 60.0;  // repeating per minute gives a quick way to assess the scan length
  redGreenBlueSpectrum(values, gradient, colour_repeat_period, replace_alpha);
}

/// write a C++ data type straight to binary format
template <typename T>
void writePlainOldData(std::ofstream &out, const T &t)
{
  out.write(reinterpret_cast<const char *>(&t), sizeof(T));
}

/// read directly from binary into a C++ data type
template <typename T>
void readPlainOldData(std::ifstream &in, T &t)
{
  in.read(reinterpret_cast<char *>(&t), sizeof(T));
}

/// write a vector of data directly to a binary file
template <typename T>
void writePlainOldDataArray(std::ofstream &out, const std::vector<T> &array)
{
  unsigned int size = (unsigned int)array.size();
  out.write(reinterpret_cast<char *>(&size), sizeof(unsigned int));
  for (unsigned int i = 0; i < size; i++) writePlainOldData(out, array[i]);
}

/// read a vector of data directly from a binary file
template <typename T>
void readPlainOldDataArray(std::ifstream &in, std::vector<T> &array)
{
  unsigned int size;
  in.read(reinterpret_cast<char *>(&size), sizeof(unsigned int));
  array.resize(size);
  for (unsigned int i = 0; i < size; i++) readPlainOldData(in, array[i]);
}

/// Log a @c std::chrono::clock::duration to an output stream.
///
/// The resulting string displays in the smallest possible unit to show three three
/// decimal places with display units ranging from seconds to nanoseconds. The table below
/// shows some example times.
///
/// Time(s)     | Display
/// ----------- | --------
/// 0.000000018 | 18ns
/// 0.000029318 | 29.318us
/// 0.0295939   | 29.593ms
/// 0.93        | 930ms
/// 15.023      | 15.023s
/// 15.000025   | 15.000s
///
/// Note that times are truncated, not rounded.
///
/// @tparam D The duration type of the form @c std::chrono::clock::duration.
/// @param out The output stream to log to.
/// @param duration The duration to convert to string.
template <typename D>
inline std::ostream &logDuration(std::ostream &out, const D &duration)
{
  const bool negative = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() < 0;
  const char *sign = (!negative) ? "" : "-";
  D abs_duration = (!negative) ? duration : duration * -1;
  auto s = std::chrono::duration_cast<std::chrono::seconds>(abs_duration).count();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(abs_duration).count();
  ms = ms % 1000;

  if (s)
  {
    out << sign << s << "." << ms << "s";
  }
  else
  {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(abs_duration).count();
    us = us % 1000;

    if (ms)
    {
      out << sign << ms << "." << us << "ms";
    }
    else
    {
      auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(abs_duration).count();
      ns = ns % 1000;

      if (us)
      {
        out << sign << us << "." << ns << "us";
      }
      else
      {
        out << sign << ns << "ns";
      }
    }
  }
  return out;
}

inline int sign(double x)
{
  return (x > 0.0) - (x < 0.0);
}

// for similar appraoch see: https://github.com/StrandedKitty/tiles-intersect/blob/master/src/index.js
template<class T> 
void walkGrid(const Eigen::Vector3d &start, const Eigen::Vector3d &end, T &object)
{
  Eigen::Vector3d direction = end - start;
  double max_length = direction.norm();
  Eigen::Vector3i p = Eigen::Vector3d(std::floor(start[0]), std::floor(start[1]), std::floor(start[2])).cast<int>();
  const Eigen::Vector3i target = Eigen::Vector3d(std::floor(end[0]), std::floor(end[1]), std::floor(end[2])).cast<int>();
  
  const Eigen::Vector3i step(sign(direction[0]), sign(direction[1]), sign(direction[2]));
  direction /= max_length;
  float eps = 1e-10; // remove tiny about so grid walking doesn't exceed its boundary
  max_length -= eps;
  Eigen::Vector3d lengths, length_delta;
  for (int j = 0; j<3; j++)
  {
    const double to = step[j] > 0 ? (double)p[j] + 1.0 - start[j] : start[j] - (double)p[j];       
    const double dir = std::max(std::numeric_limits<double>::epsilon(), std::abs(direction[j]));
    lengths[j] = to / dir;
    length_delta[j] = 1.0 / dir;
  }
  int ax = lengths[0] < lengths[1] && lengths[0] < lengths[2] ? 0 : (lengths[1] < lengths[2] ? 1 : 2);
  if (object(p, target, 0.0, lengths[ax], max_length))
  {      
    return; // only adding to one cell
  }
  
  while (lengths[ax] < max_length) 
  {
    p[ax] += step[ax];
    const double in_length = lengths[ax];
    lengths[ax] += length_delta[ax];
    ax = lengths[0] < lengths[1] && lengths[0] < lengths[2] ? 0 : (lengths[1] < lengths[2] ? 1 : 2);
    if (object(p, target, in_length, lengths[ax], max_length))
    {
      break; // only adding to one cell
    }          
  }     
}
}  // namespace ray

#endif  // RAYLIB_RAYUTILS_H
