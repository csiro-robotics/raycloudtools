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
class RAYLIB_EXPORT Argument
{
public:
  virtual bool parse(int argc, char *argv[], int &index, bool set_value) = 0; 
  virtual ~Argument() = default;
};

/// These are for fixed formats, so without - or -- prefix.
class RAYLIB_EXPORT FixedArgument : public Argument 
{
public:
  virtual ~FixedArgument() = default;
};

/// Example: "distance"
class RAYLIB_EXPORT TextArgument : public FixedArgument 
{
public:
  TextArgument(const std::string &name): name_(name) {}
  virtual ~TextArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool);
  inline const std::string &name() const { return name_; }
private:
  std::string name_;
};

/// Example: "mycloud.ply"
class RAYLIB_EXPORT FileArgument : public FixedArgument 
{
public:
  virtual ~FileArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  std::string nameStub() const { return name_.substr(0, name_.length() - 4); } // we assume a 3 letter file type
  std::string nameExt() const { return name_.substr(name_.length() - 3); }
  inline const std::string &name() const { return name_; }
private:
  std::string name_;
};

/// Numerical values
class RAYLIB_EXPORT ValueArgument : public FixedArgument 
{
public:
  virtual ~ValueArgument() = default;
};

/// Example: "4.35"
class RAYLIB_EXPORT DoubleArgument : public ValueArgument 
{
public:
  DoubleArgument();
  DoubleArgument(double min_value, double max_value) : min_value_(min_value), max_value_(max_value) {}
  virtual ~DoubleArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline double value() const { return value_; }
private:
  double value_;
  double min_value_, max_value_;
};

/// Example: "10"
class RAYLIB_EXPORT IntArgument : public ValueArgument 
{
public:
  IntArgument();
  IntArgument(int min_value, int max_value) : min_value_(min_value), max_value_(max_value) {}
  virtual ~IntArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline int value() const { return value_; }
private:
  int value_;
  int min_value_, max_value_;
};

/// Example: "1.0,2,3.26"
class RAYLIB_EXPORT Vector3dArgument : public ValueArgument 
{
public:
  Vector3dArgument();
  Vector3dArgument(double min_element_value, double max_element_value) : min_value_(min_element_value), max_value_(max_element_value) {}
  virtual ~Vector3dArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const Eigen::Vector3d &value() const { return value_; }
private:
  Eigen::Vector3d value_;
  double min_value_, max_value_;
};

/// Example: "1.0,2.4,4,-6"
class RAYLIB_EXPORT Vector4dArgument : public ValueArgument 
{
public:
  Vector4dArgument();
  Vector4dArgument(double min_element_value, double max_element_value) : min_value_(min_element_value), max_value_(max_element_value) {}
  virtual ~Vector4dArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const Eigen::Vector4d &value() const { return value_; }
private:
  Eigen::Vector4d value_;
  double min_value_, max_value_;
};

/// Parses a list of file names, e.g. "cloud1.ply cloudB.ply cloud_x.ply"
class RAYLIB_EXPORT FileArgumentList : public FixedArgument 
{
public:
  FileArgumentList(int min_number) : min_number_(min_number) {}
  virtual ~FileArgumentList() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<FileArgument> &files() const { return files_; }
private:
  std::vector<FileArgument> files_;
  int min_number_;
};

/// A choice of different keys (strings), e.g. "min"/"max"/"newest"/"oldest"
class RAYLIB_EXPORT KeyChoice : public FixedArgument 
{
public:
  KeyChoice(const std::initializer_list<std::string> &keys) : keys_(keys), selected_id_(-1) {}
  virtual ~KeyChoice() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<std::string> &keys() const { return keys_; }
  inline int selectedID() const { return selected_id_; }
  inline const std::string &selectedKey() const { return selected_key_; }
private:
  std::vector<std::string> keys_;
  int selected_id_;
  std::string selected_key_;
};

/// A choice of different key-value pairs, e.g. "pos 1,2,3" / "distance 14.2" / "num_rays 120"
class RAYLIB_EXPORT KeyValueChoice : public FixedArgument 
{
public:
  KeyValueChoice(const std::initializer_list<std::string> &keys, const std::initializer_list<ValueArgument *> &values) : 
    keys_(keys), values_(values), selected_id_(-1) {}
  virtual ~KeyValueChoice() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<std::string> &keys() const { return keys_; }
  inline const std::vector<ValueArgument *> &values() const { return values_; }
  inline int selectedID() const { return selected_id_; }
  inline const std::string &selectedKey() const { return selected_key_; }
private:
  std::vector<std::string> keys_;
  std::vector<ValueArgument *> values_; 
  int selected_id_;
  std::string selected_key_;
};

/// A choice of different value-key pairs. Usually commands defined by their units, e.g. "13.4 cm" / "12 rays" / "3.5 sigmas"
class RAYLIB_EXPORT ValueKeyChoice : public FixedArgument 
{
public:
  ValueKeyChoice(const std::initializer_list<ValueArgument *> &values, 
                 const std::initializer_list<std::string> &keys) : 
                   values_(values), keys_(keys), selected_id_(-1) {}
  virtual ~ValueKeyChoice() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<std::string> &keys() const { return keys_; }
  inline const std::vector<ValueArgument *> &values() const { return values_; }
  inline int selectedID() const { return selected_id_; }
  inline const std::string &selectedKey() const { return selected_key_; }
private:
  std::vector<ValueArgument *> values_; 
  std::vector<std::string> keys_;
  int selected_id_;
  std::string selected_key_;
};

/// For optional arguments, with the - or -- prefix
class RAYLIB_EXPORT OptionalArgument : public Argument 
{
public:
  virtual ~OptionalArgument() = default;
};

/// Optional flag, e.g. "--enable_x" or "-e"
class RAYLIB_EXPORT OptionalFlagArgument : public OptionalArgument 
{
public:
  OptionalFlagArgument(const std::string &name, char character) : 
    name_(name), character_(character), is_set_(false) {}
  virtual ~OptionalFlagArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::string &name() const { return name_; }
  inline bool isSet() const { return is_set_; }
private:
  std::string name_;
  char character_;
  bool is_set_;
};

/// Optional keyvalue pair, e.g. "--power 4.1"
struct RAYLIB_EXPORT OptionalKeyValueArgument : OptionalArgument 
{
public:
  OptionalKeyValueArgument(const std::string &name, ValueArgument *value) : 
    name_(name), value_(value), is_set_(false) {}
  virtual ~OptionalKeyValueArgument() = default;
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::string &name() const { return name_; }
private:
  std::string name_;
  ValueArgument *value_;
  bool is_set_;
};

}  // namespace ray

#endif  // RAYLIB_RAYPARSE_H
