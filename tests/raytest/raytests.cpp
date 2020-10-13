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
    for (int i = 0; i<m1.rows(); i++)
    {
      EXPECT_GT(m1[i], m2[i]-eps);
      EXPECT_LT(m1[i], m2[i]+eps);
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
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.276144, -0.0760758, 0.065631, 2.42455, 2.13738, 1.28226, 17.539, 10.1994, 0.312633, 0.314821, 0.354389, 0.987362, 0.180375, 0.164437, 0.310405, 0.111705});
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
    compareMoments(cloud.getMoments(), {-0.0860881, -0.0688598, 0.562484, 0.0215294, 0.0272777, 0.499894, -0.345469, -0.191215, 0.625452, 3.00221, 2.51452, 1.64267, 17.2104, 10.1652, 0.313403, 0.766693, 0.41599, 0.978143, 0.3197, 0.22442, 0.388217, 0.146217});
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
    compareMoments(cloud.getMoments(), {0.0999941, -0.0399831, 1.53571, 3.19998, 3.00533, 1.43086, 0.124131, -0.0325141, 1.46963, 3.09061, 2.97752, 1.49306, 38.1685, 19.0383, 0.631202, 0.458633, 0.325229, 1, 0.382899, 0.333503, 0.380153, 0});
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
    compareMoments(cloud.getMoments(), {0.788477, -0.0771538, 1.5727, 2.45926, 3.68173, 1.62275, 0.836384, -0.0622552, 1.61549, 2.35391, 3.61942, 1.5494, 31.2445, 18.0393, 0.51611, 0.501146, 0.414814, 1, 0.375121, 0.365447, 0.391637, 0});
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
    EXPECT_EQ(command("raysplit room.ply pos 0,0.1,1.5"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_outside.ply"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 1.14434e-07, 1.02466e-07, 3.37892e-08, -0.77974, 1.03139, 1.57353, 3.67521, 2.64766, 0.485084, 17.3995, 10.279, 0.311066, 0.759795, 0.425206, 0.951355, 0.321609, 0.226785, 0.39073, 0.215125});
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
    compareMoments(cloud.getMoments(), {10.1565, 19.9673, 31.7537, 3.30478, 3.04304, 1.43343, 10.191, 19.9688, 31.8099, 3.1749, 3.00455, 1.39085, 31.2445, 18.0393, 0.51611, 0.501146, 0.414814, 1, 0.375121, 0.365447, 0.391637, 0});
  }

#if RAYLIB_WITH_QHULL
  /// Creates a terrain ray cloud, then wraps it from below, comparing the mesh to the expected results
  TEST(Basic, RayWrap)
  {
    EXPECT_EQ(command("raycreate terrain 1"), 0);
    EXPECT_EQ(command("raywrap terrain.ply upwards 1.0"), 0);
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("terrain_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {0.0124983, 0.0318782, -0.00500189, 7.44097, 7.34237, 1.21885});
  }  
#endif  // RAYLIB_WITH_QHULL
} // raytest
