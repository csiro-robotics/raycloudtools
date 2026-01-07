// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raymesh.h"

#include "raylaz.h"
#include "rayply.h"
#include "raycloudwriter.h"

#include <set>

namespace ray
{
// remove additional points that are not connected to the mesh
void Mesh::reduce()
{
  std::vector<Eigen::Vector3d> verts;
  std::vector<int> new_ids(vertices_.size(), -1);
  // we do this by iterating the index list, and only adding the vertices that are in these triangles
  for (auto &ind : index_list_)
  {
    for (int i = 0; i < 3; i++)
    {
      if (new_ids[ind[i]] == -1)
      {
        new_ids[ind[i]] = (int)verts.size();
        verts.push_back(vertices_[ind[i]]);
      }
      ind[i] = new_ids[ind[i]];
    }
  }
  vertices_ = verts;
}

// convert the mesh to a height field
void Mesh::toHeightField(Eigen::ArrayXXd &field, const Eigen::Vector3d &box_min, Eigen::Vector3d box_max,
                         double width, bool fill_gaps) const
{
  double top = box_max[2];
  box_max[2] = box_min[2] + 0.5 * width;  // ensure that the grid is only 1 voxel high
  // first convert the mesh to a list of triangles, with calculated normals
  if (index_list_.empty())
  {
    std::cerr << "Error: mesh is empty" << std::endl;
  }
  std::vector<Triangle> triangles(index_list_.size());
  for (int i = 0; i < (int)index_list_.size(); i++)
  {
    Triangle &tri = triangles[i];
    for (int j = 0; j < 3; j++) 
    {
      tri.corners[j] = vertices_[index_list_[i][j]];
    }
    tri.tested = false;
    tri.normal = (tri.corners[1] - tri.corners[0]).cross(tri.corners[2] - tri.corners[0]);
  }

  // put the triangles into a grid
  Grid<Triangle *> grid(box_min, box_max, width);
  for (auto &tri : triangles)
  {
    Eigen::Vector3d tri_min = (minVector(tri.corners[0], minVector(tri.corners[1], tri.corners[2])) - box_min) / width;
    Eigen::Vector3d tri_max = (maxVector(tri.corners[0], maxVector(tri.corners[1], tri.corners[2])) - box_min) / width;
    for (int x = (int)tri_min[0]; x <= (int)tri_max[0]; x++)
    {
      for (int y = (int)tri_min[1]; y <= (int)tri_max[1]; y++) 
      {
        grid.insert(x, y, 0, &tri);
      }
    }
  }
  // now look up the triangle for each pixel centre
  const double unset = std::numeric_limits<double>::lowest();
  field = Eigen::ArrayXXd::Constant(grid.dims[0], grid.dims[1], unset);
  std::cout << "dims for low: " << grid.dims.transpose() << ", rows: " << field.rows() << ", cols: " << field.cols()
            << std::endl;
  int num_heights = 0;
  for (int x = 0; x < grid.dims[0]; x++)
  {
    for (int y = 0; y < grid.dims[1]; y++)
    {
      Eigen::Vector3d pos_top = box_min + width * (Eigen::Vector3d((double)x + 0.5, (double)y + 0.5, 0));
      Eigen::Vector3d pos_base = pos_top;
      pos_top[2] = top;
      pos_base[2] = box_min[2];
      auto &tris = grid.cell(x, y, 0).data;
      // search the triangles in this cell 'bucket'
      for (auto &tri : tris)
      {
        double depth;
        if (tri->intersectsRay(pos_top, pos_base, depth))
        {
          // intersects so interpolate the height
          double height = pos_top[2] + (pos_base[2] - pos_top[2]) * depth;
          field(x, y) = height;
          num_heights++;
          break;
        }
      }
    }
  }
  if (num_heights == 0) 
  {
    std::cout << "warning mesh does not intersect any pixel centres, using nearest triangle centre heights" << std::endl;
    for (auto &tri: triangles)
    {
      Eigen::Vector3d centre = (tri.corners[0] + tri.corners[1] + tri.corners[3])/3.0;
      Eigen::Vector3d c = (centre - box_min) / width;
      int X = std::max(0, std::min((int)std::floor(c[0]), grid.dims[0]-1)); // move to nearest inside grid
      int Y = std::max(0, std::min((int)std::floor(c[0]), grid.dims[1]-1));
      field(X,Y) = centre[2];
      num_heights++;
    }
  }
  if (num_heights == 0)
  {
    std::cerr << "Error: mesh is empty. Shouldn't get here" << std::endl;
  }
  // lastly, we repeatedly fill in the gaps
  bool gaps_remain = true;
  if (fill_gaps)
  {
    while (gaps_remain)
    {
      gaps_remain = false;
      for (int x = 0; x < grid.dims[0]; x++)
      {
        for (int y = 0; y < grid.dims[1]; y++)
        {
          if (field(x, y) == unset)
          {
            double count = 0;
            double total_height = 0;
            // look at the Moore neighbourhood to obtain a mean neighbour height
            for (int i = std::max(0, x - 1); i <= std::min(x + 1, grid.dims[0] - 1); i++)
            {
              for (int j = std::max(0, y - 1); j <= std::min(y + 1, grid.dims[1] - 1); j++)
              {
                if (field(i, j) != unset)
                {
                  total_height += field(i, j);
                  count++;
                }
              }
            }
            // Note that this immediate modifier is not order/direction independant,
            // but it doesn't matter too much, as there should be very few gaps anyway
            if (count > 0)
              field(x, y) = total_height / count;
            else
              gaps_remain = true;
          }
        }
      }
    }
  }
}

bool Mesh::splitCloud(const std::string &cloud_name, double offset, const std::string &inside_name, const std::string &outside_name)
{
  // Firstly, find the average vertex normals
  std::vector<Eigen::Vector3d> normals(vertices_.size());
  for (auto &normal : normals) normal.setZero();
  for (auto &index : index_list_)
  {
    Eigen::Vector3d normal =
      (vertices_[index[1]] - vertices_[index[0]]).cross(vertices_[index[2]] - vertices_[index[0]]);
    for (int i = 0; i < 3; i++) normals[index[i]] += normal;
  }
  for (auto &normal : normals) 
  {
    normal.normalize();
  }

  // convert to separate triangles for convenience
  std::vector<Triangle> triangles(index_list_.size());
  double mx = std::numeric_limits<double>::max();
  double mn = std::numeric_limits<double>::lowest();
  Eigen::Vector3d box_min(mx, mx, mx), box_max(mn, mn, mn);
  for (int i = 0; i < (int)index_list_.size(); i++)
  {
    Triangle &tri = triangles[i];
    for (int j = 0; j < 3; j++) tri.corners[j] = vertices_[index_list_[i][j]];
    tri.tested = false;
    tri.normal = (tri.corners[1] - tri.corners[0]).cross(tri.corners[2] - tri.corners[0]).normalized();
    for (int j = 0; j < 3; j++)
    {
      box_min = minVector(box_min, tri.corners[j]);
      box_max = maxVector(box_max, tri.corners[j]);
    }
  }

  // Thirdly, put the triangles into a grid
  double voxel_width = 1.0;
  Grid<Triangle *> grid(box_min, box_max, voxel_width);
  Grid<Triangle *> expanded_triangle_grid;

  for (auto &tri : triangles)
  {
    Eigen::Vector3d tri_min =
      (minVector(tri.corners[0], minVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
    Eigen::Vector3d tri_max =
      (maxVector(tri.corners[0], maxVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
    for (int x = (int)tri_min[0]; x <= (int)tri_max[0]; x++)
    {
      for (int y = (int)tri_min[1]; y <= (int)tri_max[1]; y++)
      {
        for (int z = (int)tri_min[2]; z <= (int)tri_max[2]; z++)
        {
          grid.insert(x, y, z, &tri);
        }
      }
    }
  }
  if (offset != 0.0)
  {
    // Thirdly, put the triangles into a grid
    expanded_triangle_grid.init(box_min, box_max, voxel_width);
    for (int i = 0; i < (int)index_list_.size(); i++)
    {
      if (!(i % 100000))
        std::cout << "filling volumes " << i << "/" << index_list_.size() << std::endl;
      Triangle &tri = triangles[i];
      Eigen::Vector3d extruded_corners[3];
      for (int j = 0; j < 3; j++) 
      {
        extruded_corners[j] = tri.corners[j] + normals[index_list_[i][j]] * offset;
      }
      Eigen::Vector3d tri_min = minVector(tri.corners[0], minVector(tri.corners[1], tri.corners[2]));
      Eigen::Vector3d tri_max = maxVector(tri.corners[0], maxVector(tri.corners[1], tri.corners[2]));
      Eigen::Vector3d tri_min2 = minVector(extruded_corners[0], minVector(extruded_corners[1], extruded_corners[2]));
      Eigen::Vector3d tri_max2 = maxVector(extruded_corners[0], maxVector(extruded_corners[1], extruded_corners[2]));
      tri_min = (minVector(tri_min, tri_min2) - box_min) / voxel_width;
      tri_max = (maxVector(tri_max, tri_max2) - box_min) / voxel_width;
      for (int x = (int)tri_min[0]; x <= (int)tri_max[0]; x++)
      {
        for (int y = (int)tri_min[1]; y <= (int)tri_max[1]; y++)
        {
          for (int z = (int)tri_min[2]; z <= (int)tri_max[2]; z++) 
          {
            expanded_triangle_grid.insert(x, y, z, &tri);
          }
        }
      }
    }
  }
  // Fourthly, drop each end point downwards to decide whether it is inside or outside..
  CloudWriter in_cloud, out_cloud;
  in_cloud.begin(inside_name);
  out_cloud.begin(outside_name);

  // splitting performed per chunk
  auto write_chunk = [&in_cloud, &out_cloud, &grid, &expanded_triangle_grid, &box_min, 
                      &voxel_width, &offset](
                    std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                    std::vector<double> &times, std::vector<RGBA> &colours) 
  {
    Cloud in_chunk, out_chunk;
    #pragma omp parallel for
    for (int i = 0; i < (int)ends.size(); i++)
    {
      int intersections = 0;
      Eigen::Vector3d start = (ends[i] - box_min) / voxel_width;
      Eigen::Vector3i index(start.cast<int>());
      std::set<Triangle *> tri_set;
      for (int z = clamped(index[2], 0, grid.dims[2] - 1); z >= 0; --z)
      {
        auto &tris = grid.cell(index[0], index[1], z).data;
        for (auto &tri : tris)
        {
          const auto &ret = tri_set.insert(tri); // downside: log n lookup, upside: allows parallelisation
          if (ret.second == false) // already exists
          {
            continue;
          }
          double depth;
          if (tri->intersectsRay(ends[i], ends[i] - Eigen::Vector3d(0.0, 0.0, 1e3), depth))
          {
            intersections++;
          }
        }
      }
      bool inside_val = offset >= 0.0;
      bool is_inside = !inside_val; // start off not inside
      if ((intersections % 2) == (int)inside_val)  // inside
      {
        bool in_tri = false;
        if (offset != 0.0) // check if it is really inside...
        {
          Eigen::Vector3d pos = (ends[i] - box_min) / voxel_width;
          Eigen::Vector3i index(pos.cast<int>());
          auto &tris = expanded_triangle_grid.cell(index[0], index[1], index[2]).data;
          for (auto &tri : tris)
          {
            if (tri->distSqrToPoint(ends[i]) < offset*offset)
            {
              in_tri = true;
              break;
            };
          }
        }
        if (offset == 0.0 || !in_tri)
        {
          is_inside = inside_val;
        }
      }
      Cloud &out = is_inside ? in_chunk : out_chunk;
      #pragma omp critical
      {
        out.addRay(starts[i], ends[i], times[i], colours[i]);
      }      
    }
    in_cloud.writeChunk(in_chunk);
    out_cloud.writeChunk(out_chunk);
  };
  if (!Cloud::read(cloud_name, write_chunk))
    return false;
  in_cloud.end();
  out_cloud.end();
  return true;
}

Eigen::Array<double, 6, 1> Mesh::getMoments() const
{
  Eigen::Array3d mean(0, 0, 0);
  for (auto &v : vertices_) 
  {
    mean += v.array();
  }
  mean /= (double)vertices_.size();
  Eigen::Array3d sigma(0, 0, 0);
  for (auto &v : vertices_) 
  {
    sigma += (v.array() - mean) * (v.array() - mean);
  }
  sigma = (sigma / (double)vertices_.size()).sqrt();
  Eigen::Array<double, 6, 1> result;
  result << mean, sigma;
  return result;
}

}  // namespace ray
