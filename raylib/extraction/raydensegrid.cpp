// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Tim Devereux
#include <nabo/nabo.h>
#include "../raycuboid.h"
#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayrenderer.h"

namespace ray
{

bool generateAreaVoxels(const std::string &cloud_stub, const double vox_width)
{
  std::string cloud_name = cloud_stub + ".ply";
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }
  const Cuboid bounds = info.ends_bound;
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  Eigen::Vector3i dims =
    (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2);  // so that we have extra space to convolve
  Cuboid grid_bounds = bounds;
  grid_bounds.min_bound_ -= Eigen::Vector3d(vox_width, vox_width, vox_width);
  DensityGrid grid(grid_bounds, vox_width, dims);
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  Eigen::MatrixXd points(grid.voxels().size(), 7);
  Eigen::MatrixXd points_q(3, grid.voxels().size());
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
          double density = grid.voxels()[index].density();
          // points_q.row(c++) = grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i+0.5, (double)j+0.5, (double)k+0.5);
          double x = grid_bounds.min_bound_[0] + vox_width * (double)i+0.5;
          double y = grid_bounds.min_bound_[1] + vox_width * (double)j+0.5;
          double z = grid_bounds.min_bound_[2] + vox_width * (double)k+0.5;
          points.row(index) << x, y, z, density, grid.voxels()[index].numHits(), grid.voxels()[index].numRays(),
            vox_width;

        }
        else
        {
          continue;
        }
      }
    }
  }

  std::ofstream myfile;
  myfile.open(cloud_stub + "_voxels.asc");

  myfile << "X Y Z DENSITY NUM_HITS NUM_RAYS VOX_WIDTH\n";

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
  printf("Wrote %d voxels to %s\n", points.rows(), (cloud_stub + "_voxels.asc").c_str());
  myfile.close();
  return true;
}


}  // namespace ray