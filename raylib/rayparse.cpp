// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayparse.h"
#include <iostream>
#include <limits>
#include "rayutils.h"

namespace ray
{
std::string getFileNameStub(const std::string &name)
{
  const size_t last_dot = name.find_last_of('.');
  if (last_dot != std::string::npos)
    return name.substr(0, last_dot);
  return name;
}

std::string getFileNameExtension(const std::string &name)
{
  const size_t last_dot = name.find_last_of('.');
  if (last_dot != std::string::npos)
    return name.substr(last_dot + 1);
  return "";
}

// Process the command line according to the specified format.
// fixed_arguments are always in order and don't have a "-" prefix. optional_arguments appear in any order after the
// fixed arguments, and have a "-" or "--" prefix.
bool parseCommandLine(int argc, char *argv[], const std::vector<FixedArgument *> &fixed_arguments,
                      std::vector<OptionalArgument *> optional_arguments, bool set_values_)
{
  // if we are setting the argument values_ then first run the parsing without setting them, and then only continue (to
  // set them) if the format matches.
  if (set_values_ && !parseCommandLine(argc, argv, fixed_arguments, optional_arguments, false))
    return false;
  int c = 1;
  for (auto &l : fixed_arguments)
  {
    if (!l->parse(argc, argv, c, set_values_))
      return false;
  }
  while (c < argc)  // now process the optional_arguments. For each command line index we need to loop through all
                    // remaining optional arguments
  {
    bool found = false;
    for (size_t i = 0; i < optional_arguments.size(); i++)
    {
      if (optional_arguments[i]->parse(argc, argv, c, set_values_))
      {
        found = true;
        optional_arguments[i] = optional_arguments.back();
        optional_arguments.pop_back();  // remove this argument, so we don't look for it twice
        break;                          // an argument should only match one Argument type
      }
    }
    if (!found)
      return false;  // no optional argument matches argument c
  }
  return true;
}

bool TextArgument::parse(int argc, char *argv[], int &index, bool)
{
  if (index >= argc)
    return false;
  bool matches = std::string(argv[index]) == name_;
  if (!matches)
    return false;
  index++;
  return true;
}

bool FileArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string file = std::string(argv[index]);
  if (file.length() < 1)
    return false;
  if (file[0] == '-') // no file should start with a dash. That is reserved for flag arguments
    return false; 
  if (check_extension_)
  {
    if (file.length() <= 2)
      return false;

    // we don't check file existence, that is up to whatever uses the file.
    // but we do check that the string contains a '.' and (if set)
    // that it has a valid >0-letter file extension
    // This lets us disambiguate files_ from other arguments
    // and isn't too restrictive, we would rather users use extensions on their file names.
    size_t extension_pos = file.rfind('.');
    if (extension_pos == std::string::npos)  // no '.' in file name
    {
      return false;
    }
    const std::string ext = file.substr(extension_pos);
    if (ext.length() < 2)
    {
      return false;
    }
    bool valid_ext = ext.at(0) == '.' && std::isalpha(ext.at(1));
    for (size_t index = 2; index < ext.length(); ++index)
    {
      valid_ext = valid_ext && std::isalnum(ext.at(index));
    }
    if (!valid_ext)
    {
      return false;
    }
  }
  if (set_value)
    name_ = file;
  index++;
  return true;
}

DoubleArgument::DoubleArgument()
{
  max_value_ = std::numeric_limits<double>::max();
  min_value_ = std::numeric_limits<double>::lowest();
}
bool DoubleArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  char *endptr;
  double val = std::strtod(argv[index], &endptr);
  if (endptr != argv[index] + std::strlen(argv[index]))  // if the double is badly formed
    return false;
  index++;
  if (!set_value)
    return true;
  bool in_range = val >= min_value_ && val <= max_value_;
  if (!in_range)
  {
    index--;
    std::cout << "Please set argument " << index << " within the range: " << min_value_ << " to " << max_value_
              << std::endl;
  }
  value_ = val;
  return in_range;
}

IntArgument::IntArgument()
{
  max_value_ = std::numeric_limits<int>::max();
  min_value_ = std::numeric_limits<int>::min();
}

bool IntArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  char *endptr;
  long int val = std::strtol(argv[index], &endptr, 10);
  if (endptr != argv[index] + std::strlen(argv[index]))  // if the int is badly formed
    return false;
  index++;
  if (!set_value)
    return true;
  bool in_range = val >= (long int)min_value_ && val <= (long int)max_value_;
  if (!in_range)
  {
    index--;
    std::cout << "Please set argument " << index << " within the range: " << min_value_ << " to " << max_value_
              << std::endl;
  }
  value_ = (int)val;
  return in_range;
}

Vector2dArgument::Vector2dArgument()
{
  max_value_ = std::numeric_limits<double>::max();
  min_value_ = std::numeric_limits<double>::lowest();
}

bool Vector2dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index]);
  std::string field;
  int i = 0;
  while (std::getline(ss, field, ','))
  {
    if (i == 2)
      return false;
    char *endptr;
    const char *str = field.c_str();
    double val = std::strtod(str, &endptr);
    if (endptr != str + std::strlen(str))  // if the double is badly formed
      return false;
    if (set_value)
    {
      if (val < min_value_ || val > max_value_)
      {
        std::cout << "Please set argument " << index << " within the range: " << min_value_ << " to " << max_value_
                  << std::endl;
        return false;
      }
      value_[i] = val;
    }
    i++;
  }
  if (i != 2)
    return false;
  index++;
  return true;
}

Vector3dArgument::Vector3dArgument()
{
  max_value_ = std::numeric_limits<double>::max();
  min_value_ = std::numeric_limits<double>::lowest();
}

bool Vector3dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index]);
  std::string field;
  int i = 0;
  while (std::getline(ss, field, ','))
  {
    if (i == 3)
      return false;
    char *endptr;
    const char *str = field.c_str();
    double val = std::strtod(str, &endptr);
    if (endptr != str + std::strlen(str))  // if the double is badly formed
      return false;
    if (set_value)
    {
      if (val < min_value_ || val > max_value_)
      {
        std::cout << "Please set argument " << index << " within the range: " << min_value_ << " to " << max_value_
                  << std::endl;
        return false;
      }
      value_[i] = val;
    }
    i++;
  }
  if (i != 3)
    return false;
  index++;
  return true;
}

Vector4dArgument::Vector4dArgument()
{
  max_value_ = std::numeric_limits<double>::max();
  min_value_ = std::numeric_limits<double>::lowest();
}

bool Vector4dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index]);
  std::string field;
  int i = 0;
  while (std::getline(ss, field, ','))
  {
    if (i == 4)
      return false;
    char *endptr;
    const char *str = field.c_str();
    double val = std::strtod(str, &endptr);
    if (endptr != str + std::strlen(str))  // if the double is badly formed
      return false;
    if (set_value)
    {
      if (val < min_value_ || val > max_value_)
      {
        std::cout << "Please set argument " << index << " within the range: " << min_value_ << " to " << max_value_
                  << std::endl;
        return false;
      }
      value_[i] = val;
    }
    i++;
  }
  if (i != 4)
    return false;
  index++;
  return true;
}

bool FileArgumentList::parse(int argc, char *argv[], int &index, bool set_value)
{
  FileArgument arg(check_extension_);
  int count = 0;
  while (arg.parse(argc, argv, index, set_value))
  {
    if (set_value)
      files_.push_back(arg);
    count++;
  }
  return count >= min_number_;
}

bool KeyChoice::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index]);
  index++;
  for (size_t i = 0; i < keys_.size(); i++)
  {
    if (keys_[i] == str)
    {
      if (set_value)
      {
        selected_id_ = (int)i;
        selected_key_ = str;
      }
      return true;
    }
  }
  return false;
}

bool KeyValueChoice::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index++]);
  for (size_t i = 0; i < keys_.size(); i++)
  {
    if (keys_[i] == str)
    {
      if (set_value)
      {
        selected_id_ = (int)i;
        selected_key_ = str;
      }
      if (!values_[i]->parse(argc, argv, index, set_value))
      {
        index--;
        return false;
      }
      return true;
    }
  }
  index--;
  return false;
}

bool ValueKeyChoice::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index + 1 >= argc)
    return false;
  std::string str(argv[index + 1]);
  for (size_t i = 0; i < keys_.size(); i++)
  {
    if (keys_[i] == str)
    {
      if (set_value)
      {
        selected_id_ = (int)i;
        selected_key_ = str;
      }
      bool parsed = values_[i]->parse(argc, argv, index, set_value);
      if (parsed)
        index++;
      return parsed;
    }
  }
  return false;
}

bool OptionalFlagArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index]);
  if (str == ("--" + name_) || str == ("-" + std::string(1, character_)))
  {
    if (set_value)
      is_set_ = true;
    index++;  // for optional parameters, we only increment the argument index when it has been found
    return true;
  }
  return false;
}

bool OptionalKeyValueArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index]);
  if (str == ("--" + name_) || str == ("-" + std::string(1, character_)))
  {
    if (set_value)
      is_set_ = true;

    index++;  // for optional parameters, we only increment the argument index when it has been found
    if (!value_->parse(argc, argv, index, set_value))
    {
      index--;
      return false;
    }
    return true;
  }
  return false;
}

}  // namespace ray
