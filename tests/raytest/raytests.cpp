// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raycloud.h"
#include "raymesh.h"
#include "rayply.h"
#include "rayforeststructure.h"
#include <vector>
#include <gtest/gtest.h>
#include <cstdlib>

/// Raycloud testing framework. In each test, the statistics of the resulting clouds are compared to the statistics
/// of the cloud when it was confirmed to be operating correctly. 
namespace raytest
{
  /// Issues the specified system command, including the required prefix on non-windows systems.
  int command(const std::string &system_command)
  {
    #ifdef _WIN32
    return system(system_command);
    #else
    return system(("./" + system_command).c_str());
    #endif // _WIN32
  }

  /// Issues the command to copy a file, which is a platform dependent system command.
  int copy(const std::string &copy_command)
  {
    #ifdef _WIN32
    return system("copy " + copy_command);
    #else
    return system(("cp " + copy_command).c_str());
    #endif // _WIN32
  }

  /// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
  /// detects differing clouds, and always equal clouds, given a tolerance @c eps.
  void compareMoments(const Eigen::ArrayXd &m1, const std::vector<double> &m2, double eps = 0.1)
  {
    for (size_t i = 0; i<m2.size(); i++)
    {
      EXPECT_GT(m1[i], m2[i]-eps);
      EXPECT_LT(m1[i], m2[i]+eps);
    }
  }

  /// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
  /// detects differing clouds, and always equal clouds, given a tolerance @c eps.
  void compareMomentsPercentageError(const Eigen::ArrayXd &m1, const std::vector<double> &m2, double percentage = 5.0)
  {
    double eps = 0.01 * percentage;
    for (size_t i = 0; i<m2.size(); i++)
    {
      EXPECT_GE(m1[i], m2[i] * (1.0 - eps));
      EXPECT_LE(m1[i], m2[i] * (1.0 + eps));
    }
  }

  /// Creates two copies of the same room with a rotational difference, then aligns the first onto the second 
  TEST(Basic, RayAlign)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(copy("room.ply room2.ply"), 0);
    EXPECT_EQ(command("rayrotate room2.ply 0,0,35"), 0);
    EXPECT_EQ(command("rayalign room.ply room2.ply"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_aligned.ply"));
    compareMoments(cloud.getMoments(), {-0.0618268, -0.077552, 0.0531072, 7.58334e-08, 7.97642e-08, 1.93877e-08, -0.180532, -0.219257, 0.0654452, 2.47241, 2.08183, 1.28226, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});  }

  /// Colours a room according to the normal direction of the surfaces, comparing to the expected results
  TEST(Basic, RayColour)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raycolour room.ply normal"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_coloured.ply"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.276144, -0.0760758, 0.065631, 2.42455, 2.13738, 1.28226, 17.539, 10.1994, 0.497919, 0.496369, 0.490293, 0.987362, 0.248361, 0.203648, 0.385192, 0.111705});
  }
  
  /// Creates two rooms, with different transformations, then combines them, and compares to the expected result.
  TEST(Basic, RayCombine)
  {
    EXPECT_EQ(command("./raycreate room 1"), 0);
    EXPECT_EQ(copy("room.ply room2.ply"), 0);
    EXPECT_EQ(command("./raytranslate room2.ply 0,0,1"), 0);
    EXPECT_EQ(command("./rayrotate room2.ply 0,0,35"), 0);
    EXPECT_EQ(command("./raycombine min room.ply room2.ply 1 rays"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_combined.ply"));
    compareMoments(cloud.getMoments(), {-0.0867714, -0.0679941, 0.546619, 0.0215326, 0.0272819, 0.499969, -0.305657, -0.186353, 0.582642, 2.95777, 2.47531, 1.63323, 17.4967, 10.1789, 0.305355, 0.763356, 0.427376, 0.979005, 0.318409, 0.225661, 0.389366, 0.143369});
  }
  
  /// Creates a building with random seed 1, and compares to the expected results
  TEST(Basic, RayCreate)
  {
    EXPECT_EQ(command("raycreate building 1"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("building.ply"));
    compareMoments(cloud.getMoments(), {-3.168, 16.5472, 7.04175, 3.28539, 18.0871, 2.96654, -3.19551, 16.564, 7.30826, 4.14359, 18.2463, 3.34431, 935.715, 540.236, 0.499997, 0.500408, 0.429416, 0.998894, 0.372134, 0.372389, 0.390499, 0.0332398});
  }
  
  /// Creates a forest and decimates it, comparing to the expected result 
  TEST(Basic, RayDecimate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("raydecimate forest.ply 10 cm"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest_decimated.ply"));
    // Below does not compare the time values (or the colour values, which are based on time here)
    // because spatial decimation does not constraint which time it picks points from.
    compareMoments(cloud.getMoments(), {-0.222571, 1.08156, 1.67264, 6.00755, 5.78731, 0.508713, -0.202668, 1.09517, 2.6238, 6.0285, 5.85715, 3.22093, 69.0574, 35.2775, 0.48969, 0.498403, 0.443549, 1, 0.379062, 0.366963, 0.389535, 0});
  }

  /// Creates a room, and calls denoise using a fixed distance threshols, and compares to expected result
  TEST(Basic, RayDenoise)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raydenoise room.ply 3 cm"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_denoised.ply"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 8.67026e-08, 8.81787e-08, 2.24394e-08, -0.464107, -0.113806, 0.161496, 2.82122, 2.34281, 1.35279, 17.81, 10.2005, 0.297047, 0.758802, 0.440232, 0.975166, 0.317215, 0.226682, 0.390971, 0.155618});
  }

  /// Creates two rooms, the second is decimated and transformed, then rayrestore is called to apply this transformation to
  /// the first (high resolution) room
  TEST(Basic, RayRestore)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(copy("room.ply room2.ply"), 0);
    EXPECT_EQ(command("raydecimate room2.ply 10 cm"), 0);
    EXPECT_EQ(command("raytranslate room2_decimated.ply 1,2,3"), 0);
    EXPECT_EQ(command("rayrotate room2_decimated.ply 0,0,-50"), 0);
    EXPECT_EQ(command("rayrestore room2_decimated.ply 10 cm room.ply"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_restored.ply"));
    compareMoments(cloud.getMoments(), {2.07399, 0.575952, 3.05217, 7.85442e-08, 7.70963e-08, 1.93877e-08, 1.9391, 0.682169, 3.06563, 2.10068, 2.45642, 1.28226, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});
  }  

  /// Creates a forest and rotates it in all three axes, comparing to the expected result
  TEST(Basic, RayRotate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("rayrotate forest.ply 10,20,30"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.ply"));
    compareMoments(cloud.getMoments(), {-0.254879, 0.846076, 2.02322, 5.62648, 5.68622, 2.56306, 0.266873, 0.772016, 3.2888, 5.54595, 5.69021, 4.28394, 62.683, 36.1903, 0.514327, 0.504407, 0.413534, 1, 0.372377, 0.365965, 0.391709, 0});
  }  

  /// Creates a room and smooths this ray cloud, comparing to the expected result
  TEST(Basic, RaySmooth)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raysmooth room.ply"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_smooth.ply"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.27615, -0.0761079, 0.0656267, 2.42413, 2.13691, 1.28163, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});
  }  

  /// Creates a room, then splits it around a plane, comparing agaisnt the expected result
  TEST(Basic, RaySplit)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raysplit room.ply plane 0,0.1,1.5"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_outside.ply"));
    compareMoments(cloud.getMoments(), {-0.467731, 1.05075, 1.43662, 2.20441, 1.60162, 0.106775, -0.77974, 1.03139, 1.57353, 3.67521, 2.64766, 0.485084, 17.3995, 10.279, 0.311066, 0.759795, 0.425206, 0.951355, 0.321609, 0.226785, 0.39073, 0.215125});
  }  

  /// Creates a room and runs raytransients, comparing the identified transients ray cloud to the expected results
  TEST(Basic, RayTransients)
  {
    EXPECT_EQ(command("raycreate room 2"), 0);
    EXPECT_EQ(command("raytransients min room.ply 1 rays"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_transient.ply"));
    compareMoments(cloud.getMoments(), {-1.05406, -0.240721, -0.0629182, 5.05649e-08, 3.32941e-08, 2.54759e-08, 0.268724, -0.136746, -0.596782, 1.04798, 0.921776, 0.527205, 32.1452, 6.7491, 0.205871, 0.395641, 0.884296, 1, 0.225501, 0.296487, 0.153923, 0});
  }  

  /// Creates a forest and translates it in all three axes, comparing to the expected result
  TEST(Basic, RayTranslate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("raytranslate forest.ply 10,20,30"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.ply"));
    compareMoments(cloud.getMoments(), {9.66298, 21.3454, 31.7177, 6.0926, 5.75511, 0.56438, 9.69155, 21.3605, 33.0883, 6.10555, 5.82564, 3.20507, 62.683, 36.1903, 0.514327, 0.504407, 0.413534, 1, 0.372377, 0.365965, 0.391709, 0});
  }

#if RAYLIB_WITH_QHULL
  /// Creates a terrain ray cloud, then wraps it from below, comparing the mesh to the expected results
  TEST(Basic, RayWrap)
  {
    EXPECT_EQ(command("raycreate terrain 1"), 0);
    EXPECT_EQ(command("raywrap terrain.ply upwards 1.0"), 0);
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("terrain_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {0.0386662, -1.52168, -0.139079, 3.30621, 3.35391, 0.705937});
  }  

  /// Tests extraction of terrain and extraction of trees
  TEST(Basic, RayExtract)
  {
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.ply"), 0);
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("forest_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {-0.00147491, -0.00191917, -0.0617946, 5.77995, 5.80266, 0.0426993});
    EXPECT_EQ(command("rayextract trees forest.ply forest_mesh.ply"), 0);

    ray::ForestStructure forest;
    EXPECT_TRUE(forest.load("forest_trees.txt"));
    compareMomentsPercentageError(forest.getMoments(), {20, 22.2172, 1058.42, 1.39416, 0.112794, 1.6403, 21918, 0, 75.75});

    EXPECT_EQ(command("rayextract forest forest.ply --ground forest_mesh.ply"), 0);

    ray::ForestStructure forest2;
    EXPECT_TRUE(forest2.load("forest_forest.txt"));
    compareMoments(forest2.getMoments(), {11, 8.43829, 586.427, 1.40054, 0.200644, 0, 16697, 8.49419, 3.09917});

    EXPECT_EQ(command("rayextract trunks forest.ply"), 0);

    ray::ForestStructure forest3;
    EXPECT_TRUE(forest3.load("forest_trunks.txt"));
    compareMoments(forest3.getMoments(), {21, 20.0797, 1124.61, 1.60427, 0.135159, 0, 0, 0, 0});
  }  
#endif  // RAYLIB_WITH_QHULL
} // raytest
