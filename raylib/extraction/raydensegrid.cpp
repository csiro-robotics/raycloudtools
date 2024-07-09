#include <nabo/nabo.h>
#include "../raycuboid.h"
#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayrenderer.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>

namespace ray
{

bool generateAreaVoxels(const std::string &cloud_stub, const double vox_width, Eigen::Vector3d user_bounds_min, Eigen::Vector3d user_bounds_max)
{
  std::string cloud_name = cloud_stub + ".ply";
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }

  Eigen::Vector3d extent;
  const Cuboid bounds = info.ends_bound;

  auto grid_bounds_min = bounds.min_bound_;
  auto grid_bounds_max = bounds.max_bound_;


  auto isZeroVector = [](const Eigen::Vector3d& vec) -> bool {
    return vec.isApprox(Eigen::Vector3d::Zero());
  };

  if (!isZeroVector(user_bounds_min))
  {
    grid_bounds_min = user_bounds_min;
  } 

  if (!isZeroVector(user_bounds_max))
  {
    grid_bounds_max = user_bounds_max;
  }

  extent = grid_bounds_max - grid_bounds_min;  

  Eigen::Vector3i dims = (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2);  // so that we have extra space to convolve
  Cuboid grid_bounds;
  if (!isZeroVector(grid_bounds_min) && !isZeroVector(grid_bounds_max))
  {
    grid_bounds.min_bound_ = grid_bounds_min - Eigen::Vector3d(vox_width, vox_width, vox_width);
    grid_bounds.max_bound_ = grid_bounds_max;
  } 
  else 
  {
    const Cuboid bounds = info.ends_bound;
    grid_bounds = bounds;
    grid_bounds.min_bound_ -= Eigen::Vector3d(vox_width, vox_width, vox_width);
  }
  extent = grid_bounds_max - grid_bounds_min;

  DensityGrid grid(grid_bounds, vox_width, dims);
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  Eigen::MatrixXd points(grid.voxels().size(), 8);
  int c = 0;
  for (int k = 0; k < dims[2]; k++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      for (int i = 0; i < dims[0]; i++)
      {
        int index = grid.getIndex(Eigen::Vector3i(i, j, k));
        if (grid.voxels()[index].numHits() > 0 && grid.voxels()[index].numRays() > 0)
        {
          double x = grid_bounds.min_bound_[0] + vox_width * (i + 0.5);
          double y = grid_bounds.min_bound_[1] + vox_width * (j + 0.5);
          double z = grid_bounds.min_bound_[2] + vox_width * (k + 0.5);
          double path_length = grid.voxels()[index].pathLength();
          double density = grid.voxels()[index].density();
          double surface_area = (density * vox_width * vox_width * vox_width);
          points.row(c++) << x, y, z, path_length, density, surface_area, grid.voxels()[index].numHits(), grid.voxels()[index].numRays(), vox_width;
        }
      }
    }
  }

  points.conservativeResize(c, points.cols());  // Resize to actual number of points

  std::ofstream myfile;
  std::ostringstream filename;
  filename << cloud_stub << "_voxels_" << vox_width << ".asc";
  myfile.open(filename.str());

  if (!myfile.is_open()) 
  {
    std::cerr << "Unable to open file " << filename.str() << std::endl;
    return false;
  }

  myfile << "X Y Z PATH_LENGTH DENSITY SURFACE_AREA NUM_HITS NUM_RAYS VOX_WIDTH\n";

  // Iterate through each row and column of the points matrix
  for (int r = 0; r < points.rows(); ++r)
  {
    for (int c = 0; c < points.cols(); ++c)
    {
      // Write the current element to the file
      myfile << points(r, c);
      if (c != points.cols() - 1)
      {
        myfile << " ";
      }
    }
    // Add a newline character at the end of the row
    myfile << "\n";
  }

  std::cout << "Wrote " << points.rows() << " voxels to " << filename.str() << std::endl;
  myfile.close();
  return true;
}

}  // namespace ray
