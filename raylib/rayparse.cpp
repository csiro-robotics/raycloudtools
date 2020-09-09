// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <iostream>
#include <limits>
#include "rayutils.h"
#include "rayparse.h"

namespace ray
{
// set values true only sets values if the command line parses correctly according to list
bool parseCommandLine(int argc, char *argv[], const std::vector<FixedArgument *> &fixed_arguments, std::vector<OptionalArgument *> optional_arguments, bool set_values)
{
  if (set_values && !parseCommandLine(argc, argv, fixed_arguments, optional_arguments, false)) // check first and exit without side effects
    return false;
  int c = 1;
  for (auto &l: fixed_arguments)
  {
    if (!l->parse(argc, argv, c, set_values))
      return false;
  }
  while (c < argc)
  {
    bool found = false;
    for (size_t i = 0; i<optional_arguments.size(); i++)
    {
      if (optional_arguments[i]->parse(argc, argv, c, set_values))
      {
        found = true;
        optional_arguments[i] = optional_arguments.back();
        optional_arguments.pop_back(); // remove this argument, so we don't look for it twice
        break; // an argument should only match one Argument type
      }
    }
    if (!found)
      return false; // no optional argument matches argument c
  }
  return true;
}

bool TextArgument::parse(int argc, char *argv[], int &index, bool)
{
  if (index >= argc)
    return false;
  return std::string(argv[index++]) == name;
}

bool FileArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string file = std::string(argv[index]);
  if (file.length() <= 4)
    return false;
  // we don't check file existence, that is up to whatever uses the file.
  // but we do check that it has a 3-letter file extension. This lets us disambiguate files from other arguments
  std::string ext = file.substr(file.length() - 4);
  bool valid_ext = ext.at(0) == '.' && std::isalpha(ext.at(1)) && std::isalnum(ext.at(2)) && std::isalnum(ext.at(3));
  if (!valid_ext) 
    return false;
  if (set_value)
    name = file;
  index++;
  return true;
}

DoubleArgument::DoubleArgument()
{ 
  max_value = std::numeric_limits<double>::max();
  min_value = std::numeric_limits<double>::lowest();
}
bool DoubleArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  char *endptr;
  double val = std::strtod(argv[index], &endptr);
  if (endptr != argv[index]+std::strlen(argv[index]))
    return false;
  index++;
  if (!set_value)
    return true;
  bool in_range = val >= min_value && val <= max_value;
  if (!in_range)
    std::cout << "Please set argument " << index << " within the range: " << min_value << " to " << max_value << std::endl;
  value = val; 
  return in_range;
}

IntArgument::IntArgument() 
{ 
  max_value = std::numeric_limits<int>::max();
  min_value = std::numeric_limits<int>::min(); 
}

bool IntArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  char *endptr;
  long int val = std::strtol(argv[index], &endptr, 10);
  if (endptr != argv[index]+std::strlen(argv[index]))
    return false;
  index++;
  if (!set_value)
    return true;
  bool in_range = val >= (long int)min_value && val <= (long int)max_value;
  if (!in_range)
    std::cout << "Please set argument " << index << " within the range: " << min_value << " to " << max_value << std::endl;
  value = (int)val;
  return in_range;
}

Vector3dArgument::Vector3dArgument() 
{ 
  max_value = std::numeric_limits<double>::max();
  min_value = std::numeric_limits<double>::lowest();
}

bool Vector3dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index++]);
  std::string field;
  int i = 0;
  while (std::getline(ss, field, ','))
  {
    if (i==3)
      return false;
    char *endptr;
    const char *str = field.c_str();
    double val = std::strtod(str, &endptr);
    if (endptr != str+std::strlen(str))
      return false;
    if (set_value)
    {
      if (val < min_value || val > max_value)
      {
        std::cout << "Please set argument " << index << " within the range: " << min_value << " to " << max_value << std::endl;
        return false;
      }
      value[i] = val;
    }
    i++;
  }
  if (i != 3)
    return false;
  return true;
}

Vector4dArgument::Vector4dArgument() 
{ 
  max_value = std::numeric_limits<double>::max();
  min_value = std::numeric_limits<double>::lowest();
}

bool Vector4dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index++]);
  std::string field;
  int i = 0;
  while (std::getline(ss, field, ','))
  {
    if (i==4)
      return false;
    char *endptr;
    const char *str = field.c_str();
    double val = std::strtod(str, &endptr);
    if (endptr != str+std::strlen(str))
      return false;
    if (set_value)
    {
      if (val < min_value || val > max_value)
      {
        std::cout << "Please set argument " << index << " within the range: " << min_value << " to " << max_value << std::endl;
        return false;
      }
      value[i] = val;
    }
    i++;
  }
  if (i != 4)
    return false;
  return true;
}

bool FileArgumentList::parse(int argc, char *argv[], int &index, bool set_value)
{
  FileArgument arg;
  int count = 0;
  while (arg.parse(argc, argv, index, set_value))
  {
    if (set_value)
      files.push_back(arg);
    count++;
  }
  return count >= min_number;
} 

bool KeyChoice::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index]);
  index++;
  for (size_t i = 0; i<keys.size(); i++)
  {
    if (keys[i] == str)
    {
      if (set_value)
      {
        selected_id = (int)i;
        selected_key = str;
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
  std::string str(argv[index]);
  index++;
  for (size_t i = 0; i<keys.size(); i++)
  {
    if (keys[i] == str)
    {
      if (set_value)
      {
        selected_id = (int)i;
        selected_key = str;
      }
      return values[i]->parse(argc, argv, index, set_value);
    }
  }
  return false;
}

bool ValueKeyChoice::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index+1 >= argc)
    return false;
  std::string str(argv[index+1]);
  for (size_t i = 0; i<keys.size(); i++)
  {
    if (keys[i] == str)
    {
      if (set_value)
      {
        selected_id = (int)i;
        selected_key = str;
      }
      bool parsed = values[i]->parse(argc, argv, index, set_value);
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
  if (str == ("--" + name) || str == ("-" + std::string(1, character)))
  {
    if (set_value)
      is_set = true;
    index++; // for optional parameters, we only increment the argument index when it has been found
    return true;
  }
  return false;
}

bool OptionalKeyValueArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::string str(argv[index]);
  if (str == ("--" + name))
  {
    if (set_value) 
      is_set = true;
    index++; // for optional parameters, we only increment the argument index when it has been found
    value->parse(argc, argv, index, set_value);
    return true;
  }
  return false;
}

}  // namespace ray
