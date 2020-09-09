// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPARSE_H
#define RAYLIB_RAYPARSE_H

#include "rayutils.h"
#include <iostream>
#include <limits>

namespace ray
{
/// Parses a command line according to a given format which can include fixed arguments and then a set of optional arguments
/// Values in the passed-in lists are only set when it returns true. This allows the function to be called multiple times for different formats
/// Only make @param set_values false if you only need to know if the format matches the arguments @param argv.
bool RAYLIB_EXPORT parseCommandLine(int argc, char *argv[], const std::vector<struct FixedArgument *> &fixed_arguments, 
                      std::vector<struct OptionalArgument *> optional_arguments = std::vector<struct OptionalArgument *>(), 
                      bool set_values = true);

/// Argument structures. These are conceptually structs, they have independent data that can be accessed and modified directly, 
/// and mainly contain a just single function for parsing
struct RAYLIB_EXPORT Argument
{
  virtual bool parse(int argc, char *argv[], int &index, bool set_value) = 0; 
};

/// These are for fixed formats, so without - or -- prefix.
struct RAYLIB_EXPORT FixedArgument : Argument 
{
};

/// Example: "distance"
struct RAYLIB_EXPORT TextArgument : FixedArgument 
{
  TextArgument(const std::string &name): name(name) {}
  std::string name;
  virtual bool parse(int argc, char *argv[], int &index, bool);
};

/// Example: "mycloud.ply"
struct RAYLIB_EXPORT FileArgument : FixedArgument 
{
  std::string name;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  std::string nameStub() { return name.substr(0, name.length() - 4); } // we assume a 3 letter file type
  std::string nameExt() { return name.substr(name.length() - 3); }
};

/// Numerical values
struct RAYLIB_EXPORT ValueArgument : FixedArgument 
{
};

/// Example: "4.35"
struct RAYLIB_EXPORT DoubleArgument : ValueArgument 
{
  DoubleArgument();
  DoubleArgument(double min_value, double max_value) : min_value(min_value), max_value(max_value) {}
  double value;
  double min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// Example: "10"
struct RAYLIB_EXPORT IntArgument : ValueArgument 
{
  IntArgument();
  IntArgument(int min_value, int max_value) : min_value(min_value), max_value(max_value) {}
  int value;
  int min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// Example: "1.0,2,3.26"
struct RAYLIB_EXPORT Vector3dArgument : ValueArgument 
{
  Vector3dArgument();
  Vector3dArgument(double min_element_value, double max_element_value) : min_value(min_element_value), max_value(max_element_value) {}
  Eigen::Vector3d value;
  double min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// Example: "1.0,2.4,4,-6"
struct RAYLIB_EXPORT Vector4dArgument : ValueArgument 
{
  Vector4dArgument();
  Vector4dArgument(double min_element_value, double max_element_value) : min_value(min_element_value), max_value(max_element_value) {}
  Eigen::Vector4d value;
  double min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// Parses a list of file names, e.g. "cloud1.ply cloudB.ply cloud_x.ply"
struct RAYLIB_EXPORT FileArgumentList : FixedArgument 
{
  FileArgumentList(int min_number) : min_number(min_number) {}

  std::vector<FileArgument> files;
  int min_number;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// A choice of different keys (strings), e.g. "min"/"max"/"newest"/"oldest"
struct RAYLIB_EXPORT KeyChoice : FixedArgument 
{
  KeyChoice(const std::initializer_list<std::string> &keys) : keys(keys), selected_id(-1) {}

  std::vector<std::string> keys;
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// A choice of different key-value pairs, e.g. "pos 1,2,3" / "distance 14.2" / "num_rays 120"
struct RAYLIB_EXPORT KeyValueChoice : FixedArgument 
{
  KeyValueChoice(const std::initializer_list<std::string> &keys, const std::initializer_list<ValueArgument *> &values) : keys(keys), values(values), selected_id(-1) {}

  std::vector<std::string> keys;
  std::vector<ValueArgument *> values; 
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// A choice of different value-key pairs. Usually commands defined by their units, e.g. "13.4 cm" / "12 rays" / "3.5 sigmas"
struct RAYLIB_EXPORT ValueKeyChoice : FixedArgument 
{
  ValueKeyChoice(const std::initializer_list<ValueArgument *> &values, const std::initializer_list<std::string> &keys) : values(values), keys(keys), selected_id(-1) {}

  std::vector<ValueArgument *> values; 
  std::vector<std::string> keys;
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// For optional arguments, with the - or -- prefix
struct RAYLIB_EXPORT OptionalArgument : Argument 
{
};

/// Optional flag, e.g. "--enable_x" or "-e"
struct RAYLIB_EXPORT OptionalFlagArgument : OptionalArgument 
{
  OptionalFlagArgument(const std::string &name, char character): name(name), character(character), is_set(false) {}
  std::string name;
  char character;
  bool is_set;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

/// Optional keyvalue pair, e.g. "--power 4.1"
struct RAYLIB_EXPORT OptionalKeyValueArgument : OptionalArgument 
{
  OptionalKeyValueArgument(const std::string &name, ValueArgument *value): name(name), value(value), is_set(false) {}
  std::string name;
  ValueArgument *value;
  bool is_set;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

}  // namespace ray

#endif  // RAYLIB_RAYPARSE_H
