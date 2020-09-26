// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raycloud.h"
#include "raymesh.h"
#include "rayply.h"
#include <vector>
#include <gtest/gtest.h>
#include <stdlib.h>

namespace raytest
{
  void compareMoments(const Eigen::ArrayXd &m1, const std::vector<double> &m2)
  {
    double eps = 0.1; 
    for (int i = 0; i<m1.rows(); i++)
    {
      EXPECT_GT(m1[i], m2[i]-eps);
      EXPECT_LT(m1[i], m2[i]+eps);
    }
  }

  TEST(Basic, RayAlign)
  {
    EXPECT_TRUE(system("scripts/align.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_aligned.ply"));
    compareMoments(cloud.getMoments(), {-0.102981, 0.882468, 2.9268, 7.37868e-08, 8.80298e-08, 4.60241e-08, 0.0901989, 0.888633, 3.037, 2.15649, 2.54802, 1.36622, 17.598, 10.2334, 0.304556, 0.760121, 0.431401, 0.982679, 0.318377, 0.227388, 0.390657, 0.130466});
  }

  TEST(Basic, RayColour)
  {
    EXPECT_TRUE(system("scripts/colour.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_coloured.ply"));
    compareMoments(cloud.getMoments(), {-0.535514, -1.21112, -0.0737656, 5.55125e-08, 1.00559e-07, 4.60241e-08, -0.404791, -1.32544, 0.0361501, 1.6056, 2.92655, 1.3663, 17.598, 10.2334, 0.342047, 0.374928, 0.371647, 0.982679, 0.176907, 0.187435, 0.294587, 0.130466});
  }
  

  TEST(Basic, RayCombine)
  {
    EXPECT_TRUE(system("scripts/combine.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_combined.ply"));
    compareMoments(cloud.getMoments(), {-0.152444, -1.25377, 0.410206, 0.395553, 0.0440426, 0.499743, 0.0743741, -1.64351, 0.594067, 2.30498, 3.27752, 1.73886, 17.4375, 10.2337, 0.3093, 0.76179, 0.424979, 0.971373, 0.319039, 0.226942, 0.390344, 0.166756});
  }
  
  TEST(Basic, RayCreate)
  {
    EXPECT_TRUE(system("scripts/create.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("building.ply"));
    compareMoments(cloud.getMoments(), {9.79481, -7.40221, -11.8606, 7.22584, 10.1319, 3.73398, 9.83437, -7.4004, -11.7155, 7.4858, 10.3845, 3.9268, 1062.12, 613.217, 0.49769, 0.503013, 0.429473, 0.999474, 0.372269, 0.373062, 0.389571, 0.0229353});
  }
  
  TEST(Basic, RayDecimate)
  {
    EXPECT_TRUE(system("scripts/decimate.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest_decimated.ply"));
    compareMoments(cloud.getMoments(), {0.401683, -0.867276, 1.58221, 3.2464, 2.85665, 1.49124, 0.446236, -0.881056, 1.49909, 3.1213, 2.79546, 1.51553, 38.9893, 19.3819, 0.628201, 0.478369, 0.318037, 1, 0.376902, 0.329625, 0.38093, 0});
  }

  TEST(Basic, RayDenoise)
  {
    EXPECT_TRUE(system("scripts/denoise.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_denoised.ply"));
    compareMoments(cloud.getMoments(), {-0.535514, -1.21112, -0.0737656, 5.39855e-08, 1.15546e-07, 5.13628e-08, -0.516675, -1.75697, 0.0589248, 1.61796, 3.38673, 1.4951, 17.8046, 10.234, 0.29876, 0.757841, 0.439485, 0.969861, 0.317083, 0.228367, 0.39102, 0.170969});
  }

  TEST(Basic, RayRestore)
  {
    EXPECT_TRUE(system("scripts/restore.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_restored.ply"));
    compareMoments(cloud.getMoments(), {0.902887, 0.151268, 2.92623, 8.47904e-08, 7.74876e-08, 4.60241e-08, 0.899335, -0.0223589, 3.03615, 2.44081, 2.27708, 1.3663, 17.598, 10.2334, 0.304556, 0.760121, 0.431401, 0.982679, 0.318377, 0.227388, 0.390657, 0.130466});
  }  

  TEST(Basic, RayRotate)
  {
    EXPECT_TRUE(system("scripts/rotate.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.ply"));
    compareMoments(cloud.getMoments(), {1.7067, -0.772091, 1.22586, 2.83001, 3.20698, 1.61407, 1.77074, -0.769162, 1.2363, 2.74548, 3.08296, 1.4453, 31.8825, 18.4077, 0.522625, 0.504213, 0.406514, 1, 0.374147, 0.362417, 0.392027, 0});
  }  

  TEST(Basic, RaySmooth)
  {
    EXPECT_TRUE(system("scripts/smooth.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_smooth.ply"));
    compareMoments(cloud.getMoments(), {-0.535514, -1.21112, -0.0737656, 5.55125e-08, 1.00559e-07, 4.60241e-08, -0.40477, -1.32547, 0.0361189, 1.60501, 2.92617, 1.36574, 17.598, 10.2334, 0.304556, 0.760121, 0.431401, 0.982679, 0.318377, 0.227388, 0.390657, 0.130466});
  }  

  TEST(Basic, RaySplit)
  {
    EXPECT_TRUE(system("scripts/split.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_outside.ply"));
    compareMoments(cloud.getMoments(), {-0.535514, -1.21112, -0.0737656, 7.3554e-08, 3.03277e-07, 1.20043e-07, 0.366716, -2.98619, 2.45637, 2.66207, 8.81238, 2.05428, 17.9612, 10.3125, 0.297631, 0.752518, 0.445934, 0.771856, 0.3174, 0.229624, 0.393306, 0.419636});
  }  

  TEST(Basic, RayTransients)
  {
    EXPECT_TRUE(system("scripts/transients.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_transient.ply"));
    compareMoments(cloud.getMoments(), {-0.24756, -0.580658, 0.503666, 6.01129e-08, 5.33004e-08, 2.95823e-08, 0.798849, -1.37795, -0.202904, 1.3397, 0.968151, 0.682122, 31.3787, 6.58515, 0.180178, 0.426694, 0.884451, 1, 0.207364, 0.297492, 0.159786, 0});
  }  

  TEST(Basic, RayTranslate)
  {
    EXPECT_TRUE(system("scripts/translate.sh"));
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.ply"));
    compareMoments(cloud.getMoments(), {10.6377, 18.8673, 31.8226, 3.32907, 2.75426, 1.49339, 10.6881, 18.8431, 31.8557, 3.16094, 2.67279, 1.41273, 31.8825, 18.4077, 0.522625, 0.504213, 0.406514, 1, 0.374147, 0.362417, 0.392027, 0});
  }  

  TEST(Basic, RayWrap)
  {
    EXPECT_TRUE(system("scripts/wrap.sh"));
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("terrain_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {0.0143732, -0.0670558, -0.0403141, 7.42256, 7.42872, 1.12602});
  }  
} // raytest
