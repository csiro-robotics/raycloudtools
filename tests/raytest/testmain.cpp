// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <gtest/gtest.h>

int main(int argc, char **argv)
{
 // raytestutils::setApplicationPath(argv[0]);
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  return err;
}
