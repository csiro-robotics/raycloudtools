// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#pragma once
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <assert.h>
#include <Eigen/Dense>

namespace RAY
{
const double pi = M_PI; 
#define ASSERT(X) assert(X);

inline std::vector<std::string> split(const std::string &s, char delim) 
{
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;
  while (getline(ss, item, delim))
    result.push_back(item);
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
T clamped(const T &value, const T &minValue, const T &maxValue)
{
  return std::max(minValue, std::min(value, maxValue));
}

template <typename T> T sgn(T val) 
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

template<typename T>
void writePlainOldData(std::ofstream &out, const T &t)
{
  out.write(reinterpret_cast<const char*>(&t), sizeof(T));
}

template<typename T>
void readPlainOldData(std::ifstream &in, T &t)
{
  in.read(reinterpret_cast<char*>(&t), sizeof(T));
}

template<typename T>
void writePlainOldDataArray(std::ofstream &out, const std::vector<T> &array)
{
  unsigned int size = array.size();
  out.write(reinterpret_cast<char*>(&size), sizeof(unsigned int));
  for (unsigned int i = 0; i<size; i++)
    writePlainOldData(out, array[i]);
}

template<typename T>
void readPlainOldDataArray(std::ifstream &in, std::vector<T> &array)
{
  unsigned int size;
  in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
  array.resize(size);
  for (unsigned int i = 0; i<size; i++)
    readPlainOldData(in, array[i]);
}

template<typename T>
void readPlainOldDataArray(std::ifstream &in, std::vector<T> &array, int decimation)
{
  unsigned int size;
  in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
  array.resize(size/decimation);
  for (unsigned int i = 0; i<size; i++)
  {
    T dat;
    readPlainOldData(in, dat);
    if ((i%decimation)==0)
      array[i/decimation] = dat;
  }
}

/// Square a value
template<class T> 
inline T sqr(const T &val)
{
  return val * val;
}

template <class T>
T mean(const std::vector<T> &list)
{
  T result = list[0];
  for (unsigned int i = 1; i<list.size(); i++)
    result += list[i];
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
  nth_element(first, middle, last); // can specify comparator as optional 4th arg
  if (list.size() % 2) // odd
    return *middle;
  else
  {
    typename std::vector<T>::iterator middle2 = middle+1;
    nth_element(first, middle2, last);
    return (*middle + *middle2)/2.0;
  }
}

/** Returns p'th percentile value in unordered list. e.g. p=50% gives median value, p=0% gives smallest value
 */
template <class T>
T percentile(std::vector<T> list, double p)
{
  typename std::vector<T>::iterator first = list.begin();
  typename std::vector<T>::iterator last = list.end();
  int closestIndex = (int)(p*((double)list.size())/100.0);
  typename std::vector<T>::iterator percentile = first + closestIndex;
  nth_element(first, percentile, last); // can specify comparator as optional 4th arg
  return *percentile;
}

// TODO: this function/macro is possibly dangerous due to its casting and manipulation of pointers. Consider alternatives.
// Retrieve a list of a certain component of a vector, e.g. a list of the w values of a vector of quaternions:
// vector<double> wValues = components(quats, quats[0].w);  // vector<Quat> quats 
// Note: this is suboptimal in that the list is passed back by value, which requires a copy. 
#define components(_list, _component) component_list(_list, _list[0]._component)
template <class U, class T>
inline std::vector <T> component_list(const std::vector<U> &list, const T &component0)
{
  unsigned long int offset = (unsigned long int)&component0 - (unsigned long int)&list[0];
  std::vector<T> subList(list.size());
  for (unsigned int i = 0; i<list.size(); ++i)
    subList[i] = *(T*)((char *)&list[i] + offset);
  return subList;
}

}