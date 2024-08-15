#include <nabo/nabo.h>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include "../raycuboid.h"
#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayrenderer.h"

namespace ray
{

bool generateDenseVoxels(const std::string &cloud_stub, const double vox_width, Eigen::Vector3d user_bounds_min,
                         Eigen::Vector3d user_bounds_max, bool write_empty)
{
  std::string cloud_name = cloud_stub + ".ply";
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }

  if (write_empty)
  {
    std::cout << "Writing empty voxels" << std::endl;
  }
  else
  {
    std::cout << "Writing non-empty voxels" << std::endl;
  }

  Eigen::Vector3d extent;
  const Cuboid bounds = info.ends_bound;
  auto grid_bounds_min = bounds.min_bound_;
  auto grid_bounds_max = bounds.max_bound_;

  auto isZeroVector = [](const Eigen::Vector3d &vec) -> bool { return vec.isApprox(Eigen::Vector3d::Zero()); };

  if (!isZeroVector(user_bounds_min))
  {
    grid_bounds_min = user_bounds_min;
  }

  if (!isZeroVector(user_bounds_max))
  {
    grid_bounds_max = user_bounds_max;
  }

  extent = grid_bounds_max - grid_bounds_min;

  // Eigen::Vector3i dims = (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2);  // so that we have extra space
  // to convolve

  // std::cout << "bounds ext: " << extent.transpose() << std::endl;
  // std::cout << "vox_width: " << vox_width << std::endl;

  // Check for division by zero
  if (std::abs(vox_width) < 1e-10)
  {
    std::cerr << "Error: vox_width is too close to zero." << std::endl;
    return false;
  }

  // Perform the calculation step by step with bounds checking
  Eigen::Vector3d temp1 = extent / vox_width;
  // std::cout << "After division: " << temp1.transpose() << std::endl;

  Eigen::Vector3d temp2 = temp1.array().ceil();
  // std::cout << "After ceiling: " << temp2.transpose() << std::endl;

  Eigen::Vector3i dims;
  for (int i = 0; i < 3; ++i)
  {
    if (temp2[i] > std::numeric_limits<int>::max() - 2)
    {
      dims[i] = std::numeric_limits<int>::max();
    }
    else
    {
      dims[i] = static_cast<int>(temp2[i]) + 2;
    }
  }
  std::cout << "Final padded dims: " << dims.transpose() << std::endl;

  // Check if dimensions are reasonable
  const int max_reasonable_dim = 1000000;  // Adjust this value as needed
  if (dims.maxCoeff() > max_reasonable_dim)
  {
    std::cerr << "Error: Resulting dimensions are unreasonably large. "
              << "Please check your input data and voxel size." << std::endl;
    return false;
  }
  Cuboid grid_bounds;
  grid_bounds.min_bound_ = grid_bounds_min - Eigen::Vector3d(vox_width, vox_width, vox_width);
  grid_bounds.max_bound_ = grid_bounds_max;

  DensityGrid grid(grid_bounds, vox_width, dims);
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  Eigen::MatrixXd points(grid.voxels().size(), 9);  // Corrected the number of columns to 9
  int c = 0;
  for (int k = 0; k < dims[2]; k++)
  {
    for (int j = 0; j < dims[1]; j++)
    {
      for (int i = 0; i < dims[0]; i++)
      {
        int index = grid.getIndex(Eigen::Vector3i(i, j, k));

        if (write_empty)
        {
          double x = grid_bounds.min_bound_[0] + vox_width * (i + 1);
          double y = grid_bounds.min_bound_[1] + vox_width * (j + 1);
          double z = grid_bounds.min_bound_[2] + vox_width * (k + 1);
          double path_length = grid.voxels()[index].pathLength();
          double density = grid.voxels()[index].density();
          double surface_area = density * vox_width * vox_width * vox_width;
          points.row(c++) << x, y, z, path_length, density, surface_area, grid.voxels()[index].numHits(),
            grid.voxels()[index].numRays(), vox_width;
        }
        else
        {
          if (grid.voxels()[index].numHits() > 0 && grid.voxels()[index].numRays() > 0)
          {
            double x = grid_bounds.min_bound_[0] + vox_width * (i + 1);
            double y = grid_bounds.min_bound_[1] + vox_width * (j + 1);
            double z = grid_bounds.min_bound_[2] + vox_width * (k + 1);
            double path_length = grid.voxels()[index].pathLength();
            double density = grid.voxels()[index].density();
            double surface_area = density * vox_width * vox_width * vox_width;
            points.row(c++) << x, y, z, path_length, density, surface_area, grid.voxels()[index].numHits(),
              grid.voxels()[index].numRays(), vox_width;
          }
        }
      }
    }
  }

  points.conservativeResize(c, points.cols());  // Resize to actual number of points

  std::ofstream myfile;
  std::ostringstream filename;
  filename << cloud_stub << "_voxels_" << vox_width << ".txt";
  myfile.open(filename.str());

  if (!myfile.is_open())
  {
    std::cerr << "Unable to open file " << filename.str() << std::endl;
    return false;
  }

  myfile << "x y z path_length density surface_area num_hits num_rays width\n";

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
