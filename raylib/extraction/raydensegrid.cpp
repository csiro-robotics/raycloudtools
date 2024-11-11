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

#if RAYLIB_WITH_NETCDF
#include <netcdf>
#endif // RAYLIB_WITH_NETCDF

namespace ray
{
bool generateDenseVoxels(const std::string &cloud_stub, const double vox_width, Eigen::Vector3d user_bounds_min,
                         Eigen::Vector3d user_bounds_max, bool write_empty, bool write_netcdf)
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

  // Check for division by zero
  if (std::abs(vox_width) < 1e-10)
  {
    std::cerr << "Error: vox_width is too close to zero." << std::endl;
    return false;
  }

  Eigen::Vector3d temp1 = extent / vox_width;
  Eigen::Vector3d temp2 = temp1.array().ceil();

  // Use Eigen::Matrix<long, 3, 1> instead of Eigen::Vector3l
  Eigen::Matrix<long, 3, 1> dims;
  for (int i = 0; i < 3; ++i)
  {
    if (temp2[i] > std::numeric_limits<long>::max() - 2)
    {
      dims[i] = std::numeric_limits<long>::max();
    }
    else
    {
      dims[i] = static_cast<long>(temp2[i]) + 2;
    }
  }
  std::cout << "Final padded dims: " << dims.transpose() << std::endl;

  // Check if dimensions are reasonable
  const long max_reasonable_dim = 1000000000L;  // 1 billion
  if (dims.maxCoeff() > max_reasonable_dim)
  {
    std::cerr << "Error: Resulting dimensions are unreasonably large. "
              << "Please check your input data and voxel size." << std::endl;
    return false;
  }

  // Calculate total number of voxels
  long long total_voxels = static_cast<long long>(dims[0]) * dims[1] * dims[2];

  // Check if total number of voxels is too large
  const long long max_total_voxels = 10000000000LL;  // 10 billion
  if (total_voxels > max_total_voxels)
  {
    std::cerr << "Error: Total number of voxels (" << total_voxels << ") exceeds maximum allowed (" << max_total_voxels
              << "). Please use a larger voxel size or smaller bounds." << std::endl;
    return false;
  }

  std::cout << "Creating Grid " << std::endl;

  Cuboid grid_bounds;
  grid_bounds.min_bound_ = grid_bounds_min - Eigen::Vector3d(vox_width, vox_width, vox_width);
  grid_bounds.max_bound_ = grid_bounds_max;
  DensityGrid grid(grid_bounds, vox_width, dims.cast<int>());

  std::cout << "Calculating Density " << std::endl;
  grid.calculateDensities(cloud_name);

  std::cout << "Adding Neighbour Priors " << std::endl;
  grid.addNeighbourPriors();

  std::cout << "Writing Grid " << std::endl;


  Eigen::MatrixXd points(grid.voxels().size(), 9);
  long c = 0;
  for (long k = 0; k < dims[2]; k++)
  {
    for (long j = 0; j < dims[1]; j++)
    {
      for (long i = 0; i < dims[0]; i++)
      {
        int index = grid.getIndex(Eigen::Vector3i(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)));

        if (write_empty)
        {
          double x = grid_bounds.min_bound_[0] + vox_width * (i + 1.5);
          double y = grid_bounds.min_bound_[1] + vox_width * (j + 1.5);
          double z = grid_bounds.min_bound_[2] + vox_width * (k + 1.5);
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
            double x = grid_bounds.min_bound_[0] + vox_width * (i + 1.5);
            double y = grid_bounds.min_bound_[1] + vox_width * (j + 1.5);
            double z = grid_bounds.min_bound_[2] + vox_width * (k + 1.5);
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

  points.conservativeResize(c, points.cols());

    if (write_netcdf) {
#if RAYLIB_WITH_NETCDF
        try {
            // Create a new NetCDF file
            std::string filename = 
                (std::ostringstream{} << cloud_stub << "_voxels_" << 
                 std::fixed << std::setprecision(1) << vox_width << ".nc").str();
            
            netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);

            // Define dimensions
            auto nPoints = dataFile.addDim("nPoints", points.rows());
            auto nVars = dataFile.addDim("nVars", points.cols());

            // Define variables
            std::vector<std::string> varNames = {
                "x", "y", "z", "path_length", "density", 
                "surface_area", "num_hits", "num_rays", "width"
            };
            
            std::vector<netCDF::NcVar> vars;
            for (const auto& varName : varNames) {
                vars.push_back(dataFile.addVar(varName, netCDF::ncDouble, {nPoints}));
            }

            // Write data
            for (size_t i = 0; i < vars.size(); ++i) {
                std::vector<double> data(points.rows());
                for (int j = 0; j < points.rows(); ++j) {
                    data[j] = points(j, i);
                }
                vars[i].putVar(data.data());
            }

            std::cout << "Wrote " << points.rows() << " voxels to " << filename << std::endl;
            return true;

        } catch (const netCDF::exceptions::NcException& e) {
            std::cerr << "NetCDF exception: " << e.what() << std::endl;
            return false;
        }
#else
        std::cerr << "NetCDF support is not enabled in this build." << std::endl;
        return false;
#endif // RAYLIB_WITH_NETCDF
    } else {
        std::ofstream myfile;
        std::ostringstream filename;
        filename << cloud_stub << "_voxels_" << vox_width << ".txt";
        
        myfile.open(filename.str());
        if (!myfile.is_open()) {
            std::cerr << "Unable to open file " << filename.str() << std::endl;
            return false;
        }

        // Write header
        myfile << "x y z path_length density surface_area num_hits num_rays width\n";

        // Write data
        for (int r = 0; r < points.rows(); ++r) {
            for (int c = 0; c < points.cols(); ++c) {
                myfile << points(r, c);
                if (c != points.cols() - 1) {
                    myfile << " ";
                }
            }
            myfile << "\n";
        }

        std::cout << "Wrote " << points.rows() << " voxels to " << filename.str() << std::endl;
        myfile.close();
        return true;
    }
}

} // namespace ray

