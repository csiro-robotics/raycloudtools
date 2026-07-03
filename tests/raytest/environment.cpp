// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Hines <thomas.hines@csiro.au>

#include "environment.hpp"

#include <filesystem>
#include <iostream>
#include <mutex>
#include <string>

namespace raytest
{
std::mutex Environment::mutex_{};
std::string Environment::path_{};

Environment::Environment(const int argc, const char * const * const argv)
{
  std::lock_guard<std::mutex> lock(mutex_);

  // If any arguments are provided, assume they are paths to the root of the
  // build directory and use them to populate the PATH environment variable so
  // that tests can execute raycloudtools executables.

  // PATH
  std::set<std::filesystem::path> path_items;
  for (int index = 0; index < argc; ++index) {
    const std::filesystem::path build_root(argv[index]);
    const auto executables_root = build_root / "raycloudtools";
    if (std::filesystem::is_directory(executables_root)) {
      for (const auto & entry : std::filesystem::directory_iterator(executables_root)) {
        if (entry.is_directory()) {
          path_items.insert(entry.path());
        }
      }
    }
  }
  if (path_items.empty()) {
    path_ = "";
  } else {
    path_ = "PATH=";
    for (const auto & path_item : path_items) {
      path_ += (path_item.string() + ":");
    }
    const auto old_path = std::getenv("PATH");
    if (old_path != nullptr) {
      path_ += std::string(old_path);
    }
    path_ += " ";
  }
}

int Environment::Command(const std::string & system_command)
{
  std::lock_guard<std::mutex> lock(mutex_);

  #ifdef _WIN32
    return system(system_command);
  #else
    const auto full_command = path_ + system_command;
    return system(full_command.c_str());
  #endif // _WIN32
}

}  // namespace raytest
