// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raycloud.h"
#include <vector>
#include <gtest/gtest.h>
#include <stdlib.h>

namespace raytest
{
  void compareMoments(const Eigen::Array<double,22,1> &m1, const std::vector<double> &m2)
  {
    double eps = 0.1; 
    for (int i = 0; i<22; i++)
    {
      EXPECT_GT(m1[i], m2[i]-eps);
      EXPECT_LT(m1[i], m2[i]+eps);
    }
  }

  TEST(Basic, RayAlign)
  {
    EXPECT_TRUE(system("scripts/align.sh"));
    ray::Cloud cloud;
    bool loaded = cloud.load("room_aligned.ply");
    EXPECT_TRUE(loaded);
    Eigen::Array<double,22,1> moments = cloud.getMoments();
    compareMoments(moments, {-0.102981, 0.882468, 2.9268, 7.37868e-08, 8.80298e-08, 4.60241e-08, 0.0901989, 0.888633, 3.037, 2.15649, 2.54802, 1.36622, 17.598, 10.2334, 0.304556, 0.760121, 0.431401, 0.982679, 0.318377, 0.227388, 0.390657, 0.130466});
  }

  TEST(Basic, RayColour)
  {
    EXPECT_TRUE(system("scripts/colour.sh"));
    ray::Cloud cloud;
    bool loaded = cloud.load("room_coloured.ply");
    EXPECT_TRUE(loaded);
    Eigen::Array<double,22,1> moments = cloud.getMoments();
    compareMoments(moments, {-0.535514, -1.21112, -0.0737656, 5.55125e-08, 1.00559e-07, 4.60241e-08, -0.404791, -1.32544, 0.0361501, 1.6056, 2.92655, 1.3663, 17.598, 10.2334, 0.342047, 0.374928, 0.371647, 0.982679, 0.176907, 0.187435, 0.294587, 0.130466});
  }
  

  TEST(Basic, RayCombine)
  {
    EXPECT_TRUE(system("scripts/combine.sh"));
    ray::Cloud cloud;
    bool loaded = cloud.load("room_combined.ply");
    EXPECT_TRUE(loaded);
    Eigen::Array<double,22,1> moments = cloud.getMoments();
    compareMoments(moments, {-0.152444, -1.25377, 0.410206, 0.395553, 0.0440426, 0.499743, 0.0743741, -1.64351, 0.594067, 2.30498, 3.27752, 1.73886, 17.4375, 10.2337, 0.3093, 0.76179, 0.424979, 0.971373, 0.319039, 0.226942, 0.390344, 0.166756});
  }
  
  TEST(Basic, RayCreate)
  {
    EXPECT_TRUE(system("scripts/create.sh"));
    ray::Cloud cloud;
    bool loaded = cloud.load("building.ply");
    EXPECT_TRUE(loaded);
    Eigen::Array<double,22,1> moments = cloud.getMoments();
    compareMoments(moments, {9.79481, -7.40221, -11.8606, 7.22584, 10.1319, 3.73398, 9.83437, -7.4004, -11.7155, 7.4858, 10.3845, 3.9268, 1062.12, 613.217, 0.49769, 0.503013, 0.429473, 0.999474, 0.372269, 0.373062, 0.389571, 0.0229353});
  }
  
  TEST(Basic, RayDecimate)
  {
    EXPECT_TRUE(system("scripts/decimate.sh"));
    ray::Cloud cloud;
    bool loaded = cloud.load("forest_decimated.ply");
    EXPECT_TRUE(loaded);
    Eigen::Array<double,22,1> moments = cloud.getMoments();
    compareMoments(moments, {0.401683, -0.867276, 1.58221, 3.2464, 2.85665, 1.49124, 0.446236, -0.881056, 1.49909, 3.1213, 2.79546, 1.51553, 38.9893, 19.3819, 0.628201, 0.478369, 0.318037, 1, 0.376902, 0.329625, 0.38093, 0});
  }
} // raytest
