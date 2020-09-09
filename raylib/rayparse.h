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
// set values true only sets values if the command line parses correctly according to list
bool parseCommandLine(int argc, char *argv[], const std::vector<struct FixedArgument *> &fixed_arguments, 
                      std::vector<struct OptionalArgument *> optional_arguments = std::vector<struct OptionalArgument *>(), 
                      bool set_values = true);


// Argument structures. These are conceptually structs, they have independent data that can be accessed and modified directly, and mainly contain a single function for parsing
struct Argument
{
  virtual bool parse(int argc, char *argv[], int &index, bool set_value) = 0; 
};

struct FixedArgument : Argument // these are for fixed formats, so without - or -- prefix.
{
};

struct TextArgument : FixedArgument // e.g. "distance"
{
  TextArgument(const std::string &name): name(name) {}
  std::string name;
  virtual bool parse(int argc, char *argv[], int &index, bool);
};

struct FileArgument : FixedArgument // e.g. "mycloud.ply"
{
  std::string name;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  std::string nameStub() { return name.substr(0, name.length() - 4); } // we assume a 3 letter file type
  std::string nameExt() { return name.substr(name.length() - 3); }
};

struct ValueArgument : FixedArgument // numerical values
{
};

struct DoubleArgument : ValueArgument // e.g. "4.35"
{
  DoubleArgument();
  DoubleArgument(double min_value, double max_value) : min_value(min_value), max_value(max_value) {}
  double value;
  double min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct IntArgument : ValueArgument // e.g. "10"
{
  IntArgument();
  IntArgument(int min_value, int max_value) : min_value(min_value), max_value(max_value) {}
  int value;
  int min_value, max_value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct Vector3dArgument : ValueArgument // e.g. "1.0,2,3.26"
{
  Eigen::Vector3d value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct Vector4dArgument : ValueArgument // e.g. "1.0,2.3,4,6"
{
  Eigen::Vector4d value;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct FileArgumentList : FixedArgument // e.g. "cloud1.ply cloudB.ply cloud_x.ply"
{
  FileArgumentList(int min_number) : min_number(min_number) {}

  std::vector<FileArgument> files;
  int min_number;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct KeyChoice : FixedArgument // e.g. "min"/"max"/"newest"/"oldest"
{
  KeyChoice(const std::initializer_list<std::string> &keys) : keys(keys), selected_id(-1) {}

  std::vector<std::string> keys;
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct KeyValueChoice : FixedArgument // e.g. "pos 1,2,3" / "distance 14.2" / "num_rays 120"
{
  KeyValueChoice(const std::initializer_list<std::string> &keys, const std::initializer_list<ValueArgument *> &values) : keys(keys), values(values), selected_id(-1) {}

  std::vector<std::string> keys;
  std::vector<ValueArgument *> values; 
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

// defined by its units (so a value-key pair)
struct ValueKeyChoice : FixedArgument // e.g. "13.4 cm" / "12 rays" / "3.5 sigmas"
{
  ValueKeyChoice(const std::initializer_list<ValueArgument *> &values, const std::initializer_list<std::string> &keys) : values(values), keys(keys), selected_id(-1) {}

  std::vector<ValueArgument *> values; 
  std::vector<std::string> keys;
  int selected_id;
  std::string selected_key;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct OptionalArgument : Argument // for optional arguments, with the - or -- prefix
{
};

struct OptionalFlagArgument : OptionalArgument // e.g. "--enable_x or -e"
{
  OptionalFlagArgument(const std::string &name, char character): name(name), character(character), is_set(false) {}
  std::string name;
  char character;
  bool is_set;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

struct OptionalKeyValueArgument : OptionalArgument // e.g. "--power 4.1"
{
  OptionalKeyValueArgument(const std::string &name, ValueArgument *value): name(name), value(value), is_set(false) {}
  std::string name;
  ValueArgument *value;
  bool is_set;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
};

}  // namespace ray

#endif  // RAYLIB_RAYPARSE_H
