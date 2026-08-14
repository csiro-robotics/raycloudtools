// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Hines <thomas.hines@csiro.au>

#ifndef RAYCLOUDTOOLS_TESTS_RAYTEST_ENVIRONMENT_HPP
#define RAYCLOUDTOOLS_TESTS_RAYTEST_ENVIRONMENT_HPP

#include <mutex>
#include <string>

#include "gtest/gtest.h"

namespace raytest
{
class Environment : public ::testing::Environment
{
public:
  Environment(const int argc, const char *const *const argv);

  /// Issues the specified system command, including the required prefix on
  /// non-Windows systems.
  static int Command(const std::string &system_command);

private:
  ///
  static std::mutex mutex_;

  ///
  static std::string path_;
};

}  // namespace raytest

#endif  // RAYCLOUDTOOLS_TESTS_RAYTEST_ENVIRONMENT_HPP
