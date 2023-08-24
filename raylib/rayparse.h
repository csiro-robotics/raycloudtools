// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPARSE_H
#define RAYLIB_RAYPARSE_H

#include <iostream>
#include <limits>
#include "rayutils.h"

namespace ray
{
/// Helper function: get the part of the filename before the .
std::string RAYLIB_EXPORT getFileNameStub(const std::string &name);

/// Helper function: get the part of the filename after the .
std::string RAYLIB_EXPORT getFileNameExtension(const std::string &name);

/// Parses a command line according to a given format which can include fixed arguments and then a set of optional
/// arguments Values in the passed-in lists are only set when it returns true. This allows the function to be called
/// multiple times for different formats Only make @param set_values false if you only need to know if the format
/// matches the arguments @param argv.
///
/// Example for command line: rayxxx raycloud.ply 0,0,30 --debug  <-- or -d, or no argument here
/// FileArgument file;
/// Vector3d argument rot(-360, 360); <-- constrain the rotation angle on each axis between -360 and 360 degrees
/// OptionalFlag debug("debug", 'd');
/// if (!ray::parseCommandLine(argc, argv, {&file, &rot}, {&debug})) <-- result in file.name() and rot.value()
///   print_usage_and_exit();
///
/// Multi-format example: 1) rayxx raycloud.ply scale 4.0   2) rayxxx smooth raycloud.ply
/// FileArgument file;
/// TextArgument scale("scale"), smooth("smooth"); <-- look for specific text
/// DoubleArgument scale_val(0.0, 100.0);  <-- viable range of values
/// bool format1 = ray::parseCommandLine(argc, argv, {&file, &scale, &scale_val}));
/// bool format2 = ray::parseCommandLine(argc, argv, {&smooth, &file}));
/// if (!format1 && !format2)
///   print_usage_and_exit();
/// Values are set only for the parseCommandLine that returned true. e.g. scale_val.value() is used if format1
bool RAYLIB_EXPORT
  parseCommandLine(int argc, char *argv[], const std::vector<struct FixedArgument *> &fixed_arguments,
                   std::vector<struct OptionalArgument *> optional_arguments = std::vector<struct OptionalArgument *>(),
                   bool set_values = true);

/// Argument structures. These are conceptually structs, they have independent data that can be accessed and modified
/// directly, and mainly contain a just single function for parsing
class RAYLIB_EXPORT Argument
{
public:
  /// Base-class parse function, called by parseCommandLine.
  /// This takes the supplied list of strings and returns true if the argument at index matches what is expected
  /// in the derived argument object. The derived object's value is only set if @param set_value is enabled.
  virtual bool parse(int argc, char *argv[], int &index, bool set_value) = 0;
  virtual ~Argument() = default;
};

/// These are for fixed formats, so without - or -- prefix.
class RAYLIB_EXPORT FixedArgument : public Argument
{
};

/// Specify a fixed piece of text, for example "range" in "raydenoise raycloud.ply range 4 cm"
class RAYLIB_EXPORT TextArgument : public FixedArgument
{
public:
  TextArgument(const std::string &name)
    : name_(name)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool);
  inline const std::string &name() const { return name_; }

private:
  std::string name_;
};

/// This is a file name (which may contain the path), it is checked that the text has an extension,
/// but the existence of the file is not checked and must be done so on any later load function
class RAYLIB_EXPORT FileArgument : public FixedArgument
{
public:
  /// @c check_extension determines whether a file's extension is checked (3 letters and alphanumeric)
  /// False is used for example for auto-merging of temporary files, which don't have standard extensions.
  FileArgument(bool check_extension = true)
    : check_extension_(check_extension)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  /// Stub is the part of the file before the '.'
  std::string nameStub() const { return getFileNameStub(name_); }
  /// The extension (excluding the dot)
  std::string nameExt() const { return getFileNameExtension(name_); }

  inline const std::string &name() const { return name_; }
  inline std::string &name() { return name_; }

private:
  std::string name_;
  bool check_extension_;
};

/// Numerical values
class RAYLIB_EXPORT ValueArgument : public FixedArgument
{
};

/// For real values, example: "4.35"
class RAYLIB_EXPORT DoubleArgument : public ValueArgument
{
public:
  DoubleArgument();
  DoubleArgument(double min_value, double max_value, double default_value = 0.0)
    : min_value_(min_value)
    , max_value_(max_value)
    , value_(default_value)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline double value() const { return value_; }

private:
  double min_value_, max_value_;
  double value_;
};

/// For integer values, example: "10"
class RAYLIB_EXPORT IntArgument : public ValueArgument
{
public:
  IntArgument();
  IntArgument(int min_value, int max_value, int default_value = 0)
    : min_value_(min_value)
    , max_value_(max_value)
    , value_(default_value)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline int value() const { return value_; }

private:
  int min_value_, max_value_;
  int value_;
};

/// For 2-component vector values, example: "1.0,2"
class RAYLIB_EXPORT Vector2dArgument : public ValueArgument
{
public:
  Vector2dArgument();
  Vector2dArgument(double min_element_value, double max_element_value, 
    const Eigen::Vector2d &default_value = Eigen::Vector2d(0,0))
    : min_value_(min_element_value)
    , max_value_(max_element_value)
    , value_(default_value)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const Eigen::Vector2d &value() const { return value_; }

private:
  double min_value_, max_value_;
  Eigen::Vector2d value_;
};

/// For 3-component vector values, example: "1.0,2,3.26"
class RAYLIB_EXPORT Vector3dArgument : public ValueArgument
{
public:
  Vector3dArgument();
  Vector3dArgument(double min_element_value, double max_element_value, 
    const Eigen::Vector3d &default_value = Eigen::Vector3d(0,0,0))
    : min_value_(min_element_value)
    , max_value_(max_element_value)
    , value_(default_value)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const Eigen::Vector3d &value() const { return value_; }

private:
  double min_value_, max_value_;
  Eigen::Vector3d value_;
};

/// For 4-component vector values, example: "1.0,2.4,4,-6"
class RAYLIB_EXPORT Vector4dArgument : public ValueArgument
{
public:
  Vector4dArgument();
  Vector4dArgument(double min_element_value, double max_element_value,
    const Eigen::Vector4d &default_value = Eigen::Vector4d(0,0,0,0))
    : min_value_(min_element_value)
    , max_value_(max_element_value)
    , value_(default_value)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const Eigen::Vector4d &value() const { return value_; }

private:
  double min_value_, max_value_;
  Eigen::Vector4d value_;
};

/// Parses a list of file names, e.g. "cloud1.ply cloudB.ply cloud_x.ply"
class RAYLIB_EXPORT FileArgumentList : public FixedArgument
{
public:
  FileArgumentList(int min_number, bool check_extension = true)
    : min_number_(min_number), check_extension_(check_extension)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<FileArgument> &files() const { return files_; }

private:
  std::vector<FileArgument> files_;
  int min_number_;
  bool check_extension_;
};

/// A choice of different keys (strings), e.g. "min"/"max"/"newest"/"oldest"
class RAYLIB_EXPORT KeyChoice : public FixedArgument
{
public:
  KeyChoice(const std::initializer_list<std::string> &keys)
    : keys_(keys)
    , selected_id_(-1)
  {}
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
  KeyValueChoice(const std::initializer_list<std::string> &keys, const std::initializer_list<FixedArgument *> &values)
    : keys_(keys)
    , values_(values)
    , selected_id_(-1)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::vector<std::string> &keys() const { return keys_; }
  inline const std::vector<FixedArgument *> &values() const { return values_; }
  inline int selectedID() const { return selected_id_; }
  inline const std::string &selectedKey() const { return selected_key_; }

private:
  std::vector<std::string> keys_;
  std::vector<FixedArgument *> values_;
  int selected_id_;
  std::string selected_key_;
};

/// A choice of different value-key pairs. Usually commands defined by their units, e.g. "13.4 cm" / "12 rays" / "3.5
/// sigmas"
class RAYLIB_EXPORT ValueKeyChoice : public FixedArgument
{
public:
  ValueKeyChoice(const std::initializer_list<ValueArgument *> &values, const std::initializer_list<std::string> &keys)
    : values_(values)
    , keys_(keys)
    , selected_id_(-1)
  {}
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
};

/// Optional flag, e.g. "--enable_x" or "-e"
class RAYLIB_EXPORT OptionalFlagArgument : public OptionalArgument
{
public:
  OptionalFlagArgument(const std::string &name, char character)
    : name_(name)
    , character_(character)
    , is_set_(false)
  {}
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
  OptionalKeyValueArgument(const std::string &name, char character, FixedArgument *value)
    : name_(name)
    , character_(character)
    , value_(value)
    , is_set_(false)
  {}
  virtual bool parse(int argc, char *argv[], int &index, bool set_value);
  inline const std::string &name() const { return name_; }
  inline bool isSet() const { return is_set_; }

private:
  std::string name_;
  char character_;
  FixedArgument *value_;
  bool is_set_;
};

}  // namespace ray

#endif  // RAYLIB_RAYPARSE_H
