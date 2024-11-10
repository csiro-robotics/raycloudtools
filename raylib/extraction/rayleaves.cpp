// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayleaves.h"
#include <nabo/nabo.h>
#include "../raycuboid.h"
#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayrenderer.h"
#define STB_IMAGE_IMPLEMENTATION
#include <random>
#include <cmath>
#include <functional>
#include <string>
#include "raylib/imageread.h"

namespace ray
{
typedef std::complex<float> Cmp;


// Function to get the leaf angle distribution based on user input
std::function<double(double)> getLeafAngleDistribution(int distribution)
{
  switch (distribution)
  {
  case 1: // Uniform distribution
    return [](double) { return 1.0; };
  case 2: // Spherical distribution
    return [](double angle) { return std::sin(angle * M_PI / 180.0); };
  case 3: // Erectophile distribution
    return [](double angle) { return std::sin(2 * angle * M_PI / 180.0); };
  case 4: // Plagiophile distribution
    return [](double angle) { return std::sin(4 * angle * M_PI / 180.0); };
  case 5: // Planophile distribution
    return [](double angle) { return std::cos(angle * M_PI / 180.0); };
  case 6: // Extremophile distribution
    return [](double angle) { return std::abs(std::cos(angle * M_PI / 180.0)); };
  default:
    return [](double) { return 1.0; };  // Default to uniform distribution
  }
}

bool generateLeaves(const std::string &cloud_stub, const std::string &trees_file, const std::string &leaf_file,
                    double leaf_area, double droop, int distribution, double leafAreaDensity, bool stalks)
{
  // For now we assume that woody points have been set as unbounded (alpha=0). e.g. through raycolour foliage or
  // raysplit file distance 0.2 as examples. so firstly we must calculate the foliage density across the whole map.
  std::string cloud_name = cloud_stub + ".ply";
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }
  const Cuboid bounds = info.ends_bound;
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  const double vox_width = 1;
  Eigen::Vector3i dims =
    (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2);  // so that we have extra space to convolve
  Cuboid grid_bounds = bounds;
  grid_bounds.min_bound_ -= Eigen::Vector3d(vox_width, vox_width, vox_width);
  DensityGrid grid(grid_bounds, vox_width, dims);
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  // we want to find the few branches that are nearest to each voxel
  // possibly we want there to be no maximum distance... which is weird, but more robust I guess.
  // so the best option is to use knn, and match voxel centres to tree segment centres I guess. This has the advantage
  // that it tends not to align leaves to really thick trunks.
  std::vector<int> tree_ids;
  std::vector<int> segment_ids;
  std::vector<std::vector<int>> neighbour_segments;  // this looks up into the above two structures
  ForestStructure forest;
  std::vector<int> dense_voxel_indices(grid.voxels().size(), -1);
  {  // Tim: this block looks for the closest cylindrical branch segments to each voxel, in order to give the leaves a
     // 'direction' value
    // The reason I use knn (K-nearest neighbour search) is that there is no maximum distance to worry about, and it is
    // fast
    if (!forest.load(trees_file))
    {
      return false;
    }

    size_t num_segments = 0;
    for (auto &tree : forest.trees)
    {
      num_segments += tree.segments().size() - 1;
    }
    size_t num_dense_voxels = 0;
    int i = 0;
    for (auto &vox : grid.voxels())
    {
      if (vox.density() > 0.0)
      {
        dense_voxel_indices[i] = (int)num_dense_voxels;
        num_dense_voxels++;
      }
      i++;
    }

    const int search_size = 12;  // find the twelve nearest branch segments. For larger voxels a larger value here would be helpful
    size_t p_size = num_segments;
    size_t q_size = num_dense_voxels;
    Eigen::MatrixXd points_p(3, p_size);
    i = 0;
    // 1. get branch centre positions
    for (int tree_id = 0; tree_id < (int)forest.trees.size(); tree_id++)
    {
      auto &tree = forest.trees[tree_id];
      for (int segment_id = 0; segment_id < (int)tree.segments().size(); segment_id++)
      {
        auto &segment = tree.segments()[segment_id];
        if (segment.parent_id == -1)
        {
          continue;
        }
        points_p.col(i++) = (segment.tip + tree.segments()[segment.parent_id].tip) / 2.0;
        tree_ids.push_back(tree_id);
        segment_ids.push_back(segment_id);
      }
    }
    // 2. get
    neighbour_segments.resize(grid.voxels().size());
    Eigen::MatrixXd points_q(3, q_size);
    int c = 0;
    for (int k = 0; k < dims[2]; k++)
    {
      for (int j = 0; j < dims[1]; j++)
      {
        for (int i = 0; i < dims[0]; i++)
        {
          int index = grid.getIndex(Eigen::Vector3i(i, j, k));
          double density = grid.voxels()[index].density();
          if (density > 0.0)
          {
            points_q.col(c++) =
              grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i + 0.5, (double)j + 0.5, (double)k + 0.5);
          }
        }
      }
    }
    Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    const double max_distance = 2.0;
    nns->knn(points_q, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, max_distance);
    delete nns;

    // Convert these set of nearest neighbours into surfels
    for (int i = 0; i < (int)grid.voxels().size(); i++)
    {
      int id = dense_voxel_indices[i];
      if (id != -1)
      {
        for (int j = 0; j < search_size && indices(j, id) != Nabo::NNSearchD::InvalidIndex; j++)
        {
          neighbour_segments[i].push_back(indices(j, id));
        }
      }
    }
  }


  // the density is now stored in grid.voxels()[grid.getIndex(Eigen::Vector3i )].density().
  struct Leaf
  {
    Eigen::Vector3d centre;
    Eigen::Vector3d direction;
    Eigen::Vector3d origin;
    double grad0;
  };
  std::vector<Leaf> leaves;
  std::vector<double> leaf_counter(grid.voxels().size());
  std::srand(1);
  for (size_t i = 0; i < grid.voxels().size(); i++)
  {
    leaf_counter[i] =
      (double)(std::rand() % 10000) / 10000.0;  // a random start stops regions of low density have 0 leaves
  }

  // Get the leaf angle distribution function based on user input
  auto leafAngleDistribution = getLeafAngleDistribution(distribution);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  auto add_leaves = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha == 0)
        continue;
      int index = grid.getIndexFromPos(ends[i]);
      auto &voxel = grid.voxels()[index];
      
      double leaf_area_per_voxel_volume = leafAreaDensity; 

      if (leaf_area_per_voxel_volume <= 0.0)
      {
        continue;
      }
      double desired_leaf_area = leaf_area_per_voxel_volume * vox_width * vox_width * vox_width;
      double num_leaves_d = desired_leaf_area / leaf_area;
      double num_points = (double)voxel.numHits();
      double &count = leaf_counter[index];
      count += num_leaves_d / num_points;
      bool add_leaf = false;
      if (count >= 1.0)
      {
        add_leaf = true;
        count--;
      }

      if (add_leaf)
      {
        Leaf new_leaf;
        new_leaf.centre = ends[i];

        double min_dist = 1e10;
        Eigen::Vector3d closest_point_on_branch(0, 0, 0);
        for (auto &ind : neighbour_segments[index])
        {
          auto &tree = forest.trees[tree_ids[ind]];
          Eigen::Vector3d line_closest;
          Eigen::Vector3d closest = tree.closestPointOnSegment(segment_ids[ind], ends[i], line_closest);
          double dist = (closest - ends[i]).norm();
          double radius = tree.segments()[segment_ids[ind]].radius;
          if (dist <= radius)
          {
            min_dist = 1e10;
            break;
          }
          if (dist < min_dist)
          {
            min_dist = dist;
            closest_point_on_branch = closest;
          }
        }
        if (min_dist == 1e10)
        {
          continue;
        }        // Calculate leaf direction using the user-specified leaf angle distribution
        Eigen::Vector3d branch_direction = (new_leaf.centre - closest_point_on_branch).normalized();
        
        // Generate a random angle using the distribution
        double angle;
        do {
          angle = dis(gen) * 90.0; // Random angle between 0 and 90 degrees
        } while (dis(gen) > leafAngleDistribution(angle));

        // Convert angle to radians
        double angle_rad = angle * M_PI / 180.0;

        // Create a rotation axis perpendicular to the branch direction
        Eigen::Vector3d rotation_axis = branch_direction.cross(Eigen::Vector3d::UnitZ()).normalized();
        if (rotation_axis.norm() < 1e-6) {
          rotation_axis = branch_direction.cross(Eigen::Vector3d::UnitY()).normalized();
        }

        // Create rotation matrix
        Eigen::AngleAxisd rotation(angle_rad, rotation_axis);
        
        // Apply rotation to branch direction to get leaf direction
        new_leaf.direction = rotation * branch_direction;

        // Apply droop
        Eigen::Vector3d flat = new_leaf.direction;
        flat[2] = 0.0;
        double dist = flat.norm();
        new_leaf.direction[2] -= droop * dist * dist;
        new_leaf.direction.normalize();

        new_leaf.origin = closest_point_on_branch;
        new_leaf.grad0 = std::tan(angle_rad);

        leaves.push_back(new_leaf);
      }  
    }
  };

  if (!ray::Cloud::read(cloud_name, add_leaves))
    return false;

  Mesh leaf_mesh;
  // could read it from file at this point
  auto &leaf_verts = leaf_mesh.vertices();
  auto &leaf_inds =
    leaf_mesh.indexList();  // one per triangle, gives the index into the vertices_ array for each corner
  auto &leaf_uvs = leaf_mesh.uvList();
  Eigen::Vector3d leaf_root(0, 0, 0);

  double leaf_width = std::sqrt(leaf_area / 2.0);
  if (leaf_file.empty())  // generate diamond leaf
  {
    // generate a 2-triangle leaf along y axis
    leaf_verts.push_back(
      Eigen::Vector3d(0, -leaf_width, -leaf_width * leaf_width * droop));  // should leaf droop just vertically?
    leaf_verts.push_back(Eigen::Vector3d(-leaf_width / 2.0, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(leaf_width / 2.0, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(0, leaf_width, -leaf_width * leaf_width * droop));
    leaf_inds.push_back(Eigen::Vector3i(0, 2, 1));
    leaf_inds.push_back(Eigen::Vector3i(2, 3, 1));
    leaf_root = leaf_verts[0];
  }
  else if (leaf_file.substr(leaf_file.length() - 4) == ".ply")  // load leaf(s) from .ply mesh
  {
    readPlyMesh(leaf_file, leaf_mesh);
    // work out its total area:
    double total_area = 0.0;
    for (auto &tri : leaf_inds)
    {
      Eigen::Vector3d side = (leaf_verts[tri[1]] - leaf_verts[tri[0]]).cross(leaf_verts[tri[2]] - leaf_verts[tri[0]]);
      total_area += side.norm() / 2.0;
    }
    double scale = std::sqrt(leaf_area / total_area);
    for (auto &vert : leaf_verts)
    {
      vert *= scale;
    }
    leaf_root = leaf_verts[0];
  }
  else if (leaf_file.substr(leaf_file.length() - 4) == ".png")  // generate leaf(s) from image
  {
    stbi_set_flip_vertically_on_load(1);
    int width, height, num_channels;
    unsigned char *image_data = stbi_load(leaf_file.c_str(), &width, &height, &num_channels, 0);
    if (!image_data)
    {
      std::cerr << "Error: cannot load file: " << leaf_file << std::endl;
      return false;
    }
    if (num_channels != 4)
    {
      std::cerr << "Error: png file has no alpha channel and no leaves are rectangles: " << leaf_file
                << ", num channels: " << num_channels << std::endl;
      return false;
    }
    double total_alpha = 0.0;
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        const int index = num_channels * (x + width * y);
        total_alpha += ((double)image_data[index + 3]) / 255.0;
      }
    }
    stbi_image_free(image_data);
    double image_area = (double)width * (double)height;
    total_alpha /= image_area;
    std::cout << "image file: " << leaf_file << " is " << total_alpha * 100.0 << "% opaque (leaf)" << std::endl;

    // generate a 4-triangle rectangle along y axis...
    // rescale the size so actual leaf coverage matches specified area
    double scale = std::sqrt(leaf_area / image_area) / total_alpha;
    double w = 0.5 * (double)width * scale;
    double h = 0.5 * (double)height * scale;
    leaf_verts.push_back(Eigen::Vector3d(-h, -w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(h, -w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(-h, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(h, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(-h, w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(h, w, -w * w * droop));
    leaf_inds.push_back(Eigen::Vector3i(0, 1, 2));
    leaf_inds.push_back(Eigen::Vector3i(2, 1, 3));
    leaf_inds.push_back(Eigen::Vector3i(2, 3, 4));
    leaf_inds.push_back(Eigen::Vector3i(4, 3, 5));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0, 0), Cmp(0, 1), Cmp(0.5, 0)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0.5, 0), Cmp(0, 1), Cmp(0.5, 1)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0.5, 0), Cmp(0.5, 1), Cmp(1, 0)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(1, 0), Cmp(0.5, 1), Cmp(1, 1)));
    leaf_mesh.textureName() = leaf_file;
    leaf_root = Eigen::Vector3d(0, -w, -w * w * droop);
  }
  else
  {
    std::cerr << "Error: leaf file type unsupported: " << leaf_file << std::endl;
    return false;
  }
  Mesh mesh;
  auto &verts = mesh.vertices();
  auto &inds = mesh.indexList();  // one per triangle, gives the index into the vertices_ array for each corner
  auto &uvs = mesh.uvList();
  mesh.textureName() = leaf_mesh.textureName();

  for (auto &leaf : leaves)
  {
    // 1. convert direction into a transformation matrix...
    Eigen::Matrix3d mat;
    mat.col(1) = leaf.direction;
    mat.col(0) = leaf.direction.cross(Eigen::Vector3d(0, 0, 1)).normalized();
    mat.col(2) = mat.col(0).cross(mat.col(1));

    int num_verts = (int)verts.size();
    for (auto &tri : leaf_inds)
    {
      inds.push_back(tri + Eigen::Vector3i(num_verts, num_verts, num_verts));
    }
    for (auto &uv : leaf_uvs)
    {
      uvs.push_back(uv);  // if UVs are present in the input, they are unchanged
    }
    for (auto &vert : leaf_verts)
    {
      verts.push_back(mat * vert + leaf.centre);
      mesh.colours().push_back(RGBA::leaves());
    }
    num_verts = (int)verts.size();
    if (stalks)
    {
      if (!uvs.empty())
      {
        std::cerr << "Error: multiple textures in one mesh are unsupported, so either turn off stalks or remove "
                     "uvs/texture from leaves"
                  << std::endl;
        return false;
      }
      Eigen::Vector3d start = leaf.origin;
      Eigen::Vector3d leaf_start = mat * leaf_root + leaf.centre;
      Eigen::Vector3d flat = (leaf_start - leaf.origin);
      flat[2] = 0.0;
      double length = flat.norm();
      flat /= length;
      Eigen::Vector3d side(-flat[1], flat[0], flat[2]);
      side *= leaf_width / 16.0;
      const int num_segs = 4;
      for (int i = 0; i < num_segs; i++)
      {
        double x = (double)i / (double)(num_segs - 1);
        x *= length;
        double h = leaf.grad0 * x - droop * x * x;
        Eigen::Vector3d pos = (i == num_segs - 1) ? leaf_start : start + Eigen::Vector3d(0, 0, h) + flat * x;
        verts.push_back(pos - side);
        verts.push_back(pos + side);
        mesh.colours().push_back(RGBA::treetrunk());
        mesh.colours().push_back(RGBA::treetrunk());
        if (i != num_segs - 1)
        {
          int j = 2 * i;
          inds.push_back(Eigen::Vector3i(num_verts, num_verts, num_verts) + Eigen::Vector3i(j, j + 2, j + 1));
          inds.push_back(Eigen::Vector3i(num_verts, num_verts, num_verts) + Eigen::Vector3i(j + 3, j + 1, j + 2));
        }
      }
    }
  }
  writePlyMesh(cloud_stub + "_leaves.ply", mesh);
  return true;
}
}  // namespace ray
