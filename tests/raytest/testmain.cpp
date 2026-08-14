// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "gtest/gtest.h"

#include "environment.h"

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new raytest::Environment(argc, argv));
  int err = RUN_ALL_TESTS();
  return err;
}
