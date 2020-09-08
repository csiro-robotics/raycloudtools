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

struct Argument
{
  virtual int numArgs() { return 1; }
  virtual bool parse(char *argv[], int &index) { return parse(argv[index++]); } 
  virtual bool parse(char *) { return false; } // this would imply not implemented by derived class
  virtual bool isOptional(){ return false; }
};

struct OptionalFlagArgument : Argument // e.g. "--enable_x or -e"
{
  OptionalFlagArgument(const std::string &name, char character): name(name), character(character), is_set(false) {}
  std::string name;
  char character;
  bool is_set;
  virtual bool parse(char *argv[], int &index) 
  {
    std::string str(argv[index]);
    if (str == ("--" + name) || str == ("-" + character))
    {
      is_set = true;
      index++; // for optional parameters, we only increment the argument index when it has been found
      return true;
    }
    return false;
  }
  virtual bool isOptional(){ return true; }
};

struct OptionalKeyValueArgument : Argument // e.g. "--power 4.1"
{
  OptionalKeyValueArgument(const std::string &name, const ValueArgument *value): name(name), value(value), is_set(false) {}
  std::string name;
  ValueArgument *value;
  bool is_set;
  virtual bool parse(char *argv[], int &index) 
  {
    std::string str(argv[index]);
    if (str == ("--" + name))
    {
      is_set = true;
      index++; // for optional parameters, we only increment the argument index when it has been found
      value->parse(argv, index);
      return true;
    }
    return false;
  }
  virtual bool isOptional(){ return true; }
  virtual int numArgs() { return 2; }
};

struct TextArgument : Argument // e.g. "distance"
{
  TextArgument(const std::string &name): name(name) {}
  std::string name;
  virtual bool parse(char *argv)
  {
    return std::string(argv) == name;
  }
};

struct FileArgument : Argument // e.g. "mycloud.ply"
{
  FileArgument() {}
  std::string name;
  virtual bool parse(char *argv)
  {
    name = std::string(argv);
    if (name.length() == 0)
      return false;
    // we don't check file existence, that is up to whatever uses the file.
    // but we do check that it has a 3-letter file extension. This lets us disambiguate files from other arguments
    std::string ext = name.substr(name.length() - 4);
    if (ext.at(0) != '.' || std::isalpha(ext.at(1))==false) // extension must be '.' followed by a letter. Can be numbers after that
      return false;
    return true;
  }
  std::string nameStub()
  {
    return name.substr(0, name.length() - 4); // we assume a 3 letter file type   
  }
  std::string nameExt()
  {
    return name.substr(name.length() - 3);    
  }
};

struct ValueArgument : Argument // numerical values
{
};

struct DoubleArgument : ValueArgument // e.g. "4.35"
{
  DoubleArgument()
  { 
    max_value = std::numeric_limits<double>::max();
    min_value = std::numeric_limits<double>::lowest();
  }
  DoubleArgument(double min_value, double max_value) : min_value(min_value), max_value(max_value) {}
  double value;
  double min_value, max_value;
  virtual bool parse(char *argv)
  {
    value = std::stod(argv);
    bool in_range = value >= min_value && value <= max_value;
    if (!in_range)
      std::cout << "Please set the command line value " << value << " within the range: " << min_value << " to " << max_value << std::endl;
    return in_range;
  }
};

struct IntArgument : ValueArgument // e.g. "10"
{
  IntArgument() 
  { 
    max_value = std::numeric_limits<int>::max();
    min_value = std::numeric_limits<int>::min(); 
  }
  IntArgument(int min_value, int max_value) : min_value(min_value), max_value(max_value) {}
  int value;
  int min_value, max_value;
  virtual bool parse(char *argv)
  {
    value = std::stoi(argv);
    bool in_range = value >= min_value && value <= max_value;
    if (!in_range)
      std::cout << "Please set the command line value " << value << " within the range: " << min_value << " to " << max_value << std::endl;
    return in_range;
  }
};

struct Vector3dArgument : ValueArgument // e.g. "1.0,2,3.26"
{
  Eigen::Vector3d value;
  virtual bool parse(char *argv)
  {
    std::stringstream ss(argv);
    ss >> value[0];
    ss.ignore(1);
    ss >> value[1];
    ss.ignore(1);
    ss >> value[2];
    return true;
  }
};

struct Vector4dArgument : ValueArgument // e.g. "1.0,2.3,4,6"
{
  Eigen::Vector4d value;
  virtual bool parse(char *argv)
  {
    std::stringstream ss(argv);
    ss >> value[0];
    ss.ignore(1);
    ss >> value[1];
    ss.ignore(1);
    ss >> value[2];
    ss.ignore(1);
    ss >> value[3];
    return true;
  }
};

struct FileArgumentList : Argument // e.g. "cloud1.ply cloudB.ply cloud_x.ply"
{
  FileArgumentList(int min_number) : min_number(min_number) {}

  std::vector<FileArgument> files;
  int min_number;
  virtual int numArgs() { return 100; }
  virtual bool parse(char *argv[], int &index) 
  {
    FileArgument arg;
    while (arg.parse(argv[index]))
    {
      files.push_back(arg);
      index++;
    }
    return (int)files.size() >= min_number;
  } 
};

struct KeyChoice : Argument // e.g. "min"/"max"/"newest"/"oldest"
{
  KeyChoice(const std::initializer_list<std::string> &keys) : keys(keys), selected_id(-1) {}

  std::vector<std::string> keys;
  int selected_id;
  std::string selected_key;
  virtual int numArgs() { return 1; }
  virtual bool parse(char *argv[], int &index)
  {
    std::string str(argv[index]);
    index++;
    for (size_t i = 0; i<keys.size(); i++)
    {
      if (keys[i] == str)
      {
        selected_id = (int)i;
        selected_key = str;
        return true;
      }
    }
    return false;
  }
};

struct KeyValueChoice : Argument // e.g. "pos 1,2,3" / "distance 14.2" / "num_rays 120"
{
  KeyValueChoice(const std::initializer_list<std::string> &keys, const std::initializer_list<ValueArgument *> &values) : keys(keys), values(values), selected_id(-1) {}

  std::vector<std::string> keys;
  std::vector<ValueArgument *> values; 
  int selected_id;
  std::string selected_key;
  virtual int numArgs() { return 2; }
  virtual bool parse(char *argv[], int &index)
  {
    std::string str(argv[index]);
    index++;
    for (size_t i = 0; i<keys.size(); i++)
    {
      if (keys[i] == str)
      {
        selected_id = (int)i;
        selected_key = str;
        return values[i]->parse(argv, index);
      }
    }
    return false;
  }
};

// defined by its units (so a value-key pair)
struct ValueKeyChoice : Argument // e.g. "13.4 cm" / "12 rays" / "3.5 sigmas"
{
  ValueKeyChoice(const std::initializer_list<ValueArgument *> &values, const std::initializer_list<std::string> &keys) : keys(keys), values(values), selected_id(-1) {}

  std::vector<std::string> keys;
  std::vector<ValueArgument *> values; // may be empty
  int selected_id;
  std::string selected_key;
  virtual int numArgs() { return 2; }
  virtual bool parse(char *argv[], int &index)
  {
    std::string str(argv[index+1]);
    for (size_t i = 0; i<keys.size(); i++)
    {
      if (keys[i] == str)
      {
        selected_id = (int)i;
        selected_key = str;
        bool parsed = values[i]->parse(argv, index);
        index++;
        return parsed;
      }
    }
    return false;
  }
};

bool parseCommandLine(int argc, char *argv[], const std::initializer_list<Argument *> &list)
{
  int minArgs = 1, maxArgs = 1;
  std::vector<Argument *> fixed_list;
  std::vector<Argument *> optional_list;
  for (auto &l: list)
  {
    if (l->isOptional())
      optional_list.push_back(l);
    else if (!optional_list.empty())
    {
      std::cout << "Error: parseCommandLine requirement- set optional arguments after fixed arguments." << std::endl;
      return false;
    }
    else 
    {
      fixed_list.push_back(l);
      minArgs += l->numArgs();
    }
    maxArgs += l->numArgs();
  }
  if (argc < minArgs || argc > maxArgs)
    return false;
  int c = 1;
  for (auto &l: fixed_list)
  {
    if (!l->parse(argv, c))
      return false;
  }
  while (c < argc)
  {
    bool found = false;
    for (size_t i = 0; i<optional_list.size(); i++)
    {
      if (optional_list[i]->parse(argv, c))
      {
        found = true;
        optional_list[i] = optional_list.back();
        optional_list.pop_back(); // remove this argument, so we don't look for it twice
        break; // an argument should only match one Argument type
      }
    }
    if (!found)
      return false; // no optional argument matches argument c
  }
  return true;
}


}  // namespace ray

#endif  // RAYLIB_RAYPARSE_H
