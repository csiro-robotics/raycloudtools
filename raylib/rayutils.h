// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYUTILS_H
#define RAYLIB_RAYUTILS_H

#include "raylib/raylibconfig.h"

#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <assert.h>
#include <set>

#include <Eigen/Dense>

namespace ray
{
const double kPi = M_PI;
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
  return min + (max - min) * (double(rand()) / double(RAND_MAX));
}

inline std::vector<int64_t> voxelSubsample(const std::vector<Eigen::Vector3d> &points, double voxel_width)
{
  struct Vector3iLess
  {
    bool operator()(const Eigen::Vector3i &a, const Eigen::Vector3i &b) const
    {
      if (a[0] != b[0])
        return a[0] < b[0];
      if (a[1] != b[1])
        return a[1] < b[1];
      return a[2] < b[2];
    }
  };
  std::vector<int64_t> indices;
  std::set<Eigen::Vector3i, Vector3iLess> test_set;
  for (int64_t i = (int64_t)points.size()-1; i >= 0; i--)
  {
    Eigen::Vector3i place(int(std::floor(points[i][0] / voxel_width)), int(std::floor(points[i][1] / voxel_width)),
                          int(std::floor(points[i][2] / voxel_width)));
    if (test_set.find(place) == test_set.end())
    {
      test_set.insert(place);
      indices.push_back(i);
    }
  }
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
  Eigen::Vector3d purple(0.25, -0.5, 0.25);
  Eigen::Vector3d orange(0.5, 0.0, -0.5);
  orange.normalize();
  orange *= purple.norm();

  gradient.resize(values.size());
  for (unsigned int i = 0; i < values.size(); i++)
  {
    double angle = (2.0 * kPi / 6.0) * (0.5 + 5.0 * (values[i] - min_value) / (max_value - min_value));
    Eigen::Vector3d col = Eigen::Vector3d(0.5, 0.5, 0.5) + purple * cos(angle) + orange * sin(angle);
    gradient[i].red = uint8_t(255.0 * (0.5 + 0.5 * col[0]));
    gradient[i].green = uint8_t(255.0 * (0.5 + 0.5 * col[1]));
    gradient[i].blue = uint8_t(255.0 * (0.5 + 0.5 * col[2]));
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline void redGreenBlueSpectrum(const std::vector<double> &values, std::vector<RGBA> &gradient, double wavelength,
                                 bool replace_alpha)
{
  Eigen::Vector3d purple(0.25, -0.5, 0.25);
  Eigen::Vector3d orange(0.5, 0.0, -0.5);
  orange.normalize();
  orange *= purple.norm();

  gradient.resize(values.size());
  for (unsigned int i = 0; i < values.size(); i++)
  {
    double angle = (2.0 * kPi) * values[i] / wavelength;
    Eigen::Vector3d col = Eigen::Vector3d(0.5, 0.5, 0.5) + purple * cos(angle) + orange * sin(angle);
    gradient[i].red = uint8_t(255.0 * (0.5 + 0.5 * col[0]));
    gradient[i].green = uint8_t(255.0 * (0.5 + 0.5 * col[1]));
    gradient[i].blue = uint8_t(255.0 * (0.5 + 0.5 * col[2]));
    if (replace_alpha)
      gradient[i].alpha = 255;
  }
}

inline void colourByTime(const std::vector<double> &values, std::vector<RGBA> &gradient, bool replace_alpha = true)
{
  redGreenBlueSpectrum(values, gradient, 10.0, replace_alpha);
}
}  // namespace ray

#endif  // RAYLIB_RAYUTILS_H
