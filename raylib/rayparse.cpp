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
  std::string file = std::string(argv[index++]);
  if (file.length() == 0)
    return false;
  // we don't check file existence, that is up to whatever uses the file.
  // but we do check that it has a 3-letter file extension. This lets us disambiguate files from other arguments
  std::string ext = file.substr(file.length() - 4);
  if (ext.at(0) != '.' || std::isalpha(ext.at(1))==false) // extension must be '.' followed by a letter. Can be numbers after that
    return false;
  if (set_value)
    name = file;
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
  double val = std::stod(argv[index++]);
  value = val; 
  if (!set_value)
    return true;
  bool in_range = val >= min_value && val <= max_value;
  if (!in_range)
    std::cout << "Please set the command line value " << value << " within the range: " << min_value << " to " << max_value << std::endl;
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
  int val = std::stoi(argv[index++]);
  value = val;
  if (!set_value)
    return true;
  bool in_range = val >= min_value && val <= max_value;
  if (!in_range)
    std::cout << "Please set the command line value " << value << " within the range: " << min_value << " to " << max_value << std::endl;
  return in_range;
}

bool Vector3dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index++]);
  if (set_value)
  {
    ss >> value[0];
    ss.ignore(1);
    ss >> value[1];
    ss.ignore(1);
    ss >> value[2];
  }
  return true;
}

bool Vector4dArgument::parse(int argc, char *argv[], int &index, bool set_value)
{
  if (index >= argc)
    return false;
  std::stringstream ss(argv[index++]);
  if (set_value)
  {
    ss >> value[0];
    ss.ignore(1);
    ss >> value[1];
    ss.ignore(1);
    ss >> value[2];
    ss.ignore(1);
    ss >> value[3];
  }
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
    index++;
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
  if (str == ("--" + name) || str == ("-" + character))
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
