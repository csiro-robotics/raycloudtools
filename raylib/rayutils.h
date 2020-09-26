// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYUTILS_H
#define RAYLIB_RAYUTILS_H

#include "raylib/raylibconfig.h"

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

#include <Eigen/Dense>
#include "rayrandom.h"

namespace ray
{
const double kPi = M_PI;
// while this is an absolute value, it has little effect on results unless point spacing is signficantly less
// than this small value in metres. However, the computation time is better for having a value greater than 0..
const double kNearestNeighbourEpsilon = 0.001; 
#define ASSERT(X) assert(X);

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

inline std::vector<int64_t> voxelSubsample(const std::vector<Eigen::Vector3d> &points, double voxel_width, std::set<Eigen::Vector3i, Vector3iLess> *voxel_set = NULL)
{
  std::vector<int64_t> indices;
  std::set<Eigen::Vector3i, Vector3iLess> vox_set;
  // We do a backwards iteration of the list in order to pick the last (newest) point in each voxel.
  // The choice of which opint per voxel is ambiguous, but while we are ordering chronologically, it benefits
  // merge-and-decimate operations to pick the most recent point
  for (int64_t i = (int64_t)points.size()-1; i >= 0; i--)
  {
    Eigen::Vector3i voxel(int(std::floor(points[i][0] / voxel_width)), int(std::floor(points[i][1] / voxel_width)),
                          int(std::floor(points[i][2] / voxel_width)));
    if (vox_set.find(voxel) == vox_set.end())
    {
      vox_set.insert(voxel);
      indices.push_back(i);
    }
  }
  // we need to reverse the list to maintain the ordering
  int64_t numIndices = (int64_t)indices.size();
  for (int64_t i = 0; i<numIndices/2; i++)
    std::swap(indices[i], indices[numIndices - 1 - i]);
  if (voxel_set != NULL)
    *voxel_set = std::move(vox_set);
  return indices;
}

/// Square a value
template <class T>
inline T sqr(const T &val)
{
  return val * val;
}

template <class T>
T mean(const std::vector<T> &list)
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
T median(std::vector<T> list)
{
  typename std::vector<T>::iterator first = list.begin();
  typename std::vector<T>::iterator last = list.end();
  typename std::vector<T>::iterator middle = first + ((last - first) / 2);
  nth_element(first, middle, last);  // can specify comparator as optional 4th arg
  if (list.size() % 2)               // odd
    return *middle;
  else
  {
    typename std::vector<T>::iterator middle2 = middle + 1;
    nth_element(first, middle2, last);
    return (*middle + *middle2) / 2.0;
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
  uint8_t red;
  uint8_t green;
  uint8_t blue;
  uint8_t alpha;
};

inline void redGreenBlueGradient(const std::vector<double> &values, std::vector<RGBA> &gradient, double min_value,
                                 double max_value, bool replace_alpha)
{
  gradient.resize(values.size());
  Eigen::Vector3d hue_cycle[6] = {Eigen::Vector3d(1.0, 0.0, 0.5), Eigen::Vector3d(1.0, 0.5, 0.0), 
                                  Eigen::Vector3d(0.5, 1.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.5), 
                                  Eigen::Vector3d(0.0, 0.5, 1.0), Eigen::Vector3d(0.5, 0.0, 1.0)};
  for (size_t i = 0; i < values.size(); i++)
  {
    double v = 0.5 + 4.0*clamped((values[i] - min_value)/(max_value - min_value), 0.0, 1.0);
    int id = (int)v;
    double blend = v - (double)id;
    Eigen::Vector3d col = hue_cycle[id]*(1.0-blend) + hue_cycle[id+1]*blend;
    gradient[i].red = uint8_t(255.0 * col[0]);
    gradient[i].green = uint8_t(255.0 * col[1]);
    gradient[i].blue = uint8_t(255.0 * col[2]);
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline void redGreenBlueSpectrum(const std::vector<double> &values, std::vector<RGBA> &gradient, double wavelength,
                                 bool replace_alpha)
{
  gradient.resize(values.size());
  Eigen::Vector3d cycle[6] = {Eigen::Vector3d(1.0, 0.5, 0.0), Eigen::Vector3d(0.5, 1.0, 0.0), 
                              Eigen::Vector3d(0.0, 1.0, 0.5), Eigen::Vector3d(0.0, 0.5, 1.0),
                              Eigen::Vector3d(0.5, 0.0, 1.0), Eigen::Vector3d(1.0, 0.0, 0.1)};
  for (size_t i = 0; i < values.size(); i++)
  {
    double temp = values[i]/wavelength;
    double v = 6.0 * (temp - std::floor(temp));
    int id = (int)v;
    double blend = v - (double)id;
    Eigen::Vector3d col = cycle[id]*(1.0-blend) + cycle[(id+1)%6]*blend;
    gradient[i].red = uint8_t(255.0 * col[0]);
    gradient[i].green = uint8_t(255.0 * col[1]);
    gradient[i].blue = uint8_t(255.0 * col[2]);
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline void colourByTime(const std::vector<double> &values, std::vector<RGBA> &gradient, bool replace_alpha = true)
{
  const double colour_repeat_period = 60.0; // repeating per minute gives a quick way to assess the scan length
  redGreenBlueSpectrum(values, gradient, colour_repeat_period, replace_alpha);
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
    out << sign << s << "." << std::setw(3) << std::setfill('0') << ms << "s";
  }
  else
  {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(abs_duration).count();
    us = us % 1000;

    if (ms)
    {
      out << sign << ms << "." << std::setw(3) << std::setfill('0') << us << "ms";
    }
    else
    {
      auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(abs_duration).count();
      ns = ns % 1000;

      if (us)
      {
        out << sign << us << "." << std::setw(3) << std::setfill('0') << ns << "us";
      }
      else
      {
        out << sign << ns << "ns";
      }
    }
  }
  return out;
}
}  // namespace ray

#endif  // RAYLIB_RAYUTILS_H
