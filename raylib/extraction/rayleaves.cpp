// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayleaves.h"
#include "../rayrenderer.h"
#include "../raycuboid.h"
#include "../rayply.h"
#include "../raymesh.h"
#include "../rayforeststructure.h"
#include <nabo/nabo.h>

namespace ray
{

bool generateLeaves(const std::string &cloud_stub, const std::string &trees_file, const std::string &leaf_file, double leaf_area)
{
  // For now we assume that woody points have been set as unbounded (alpha=0). e.g. through raycolour foliage or raysplit file distance 0.2 as examples.
  // so firstly we must calculate the foliage density across the whole map. 
  std::string cloud_name = cloud_stub + ".ply";
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }
  const Cuboid bounds = info.ends_bound; 
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  const double vox_width = 0.5;
  Eigen::Vector3i dims = (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2); // so that we have extra space to convolve
  Cuboid grid_bounds = bounds;
  grid_bounds.min_bound_ -= Eigen::Vector3d(vox_width, vox_width, vox_width);
  DensityGrid grid(grid_bounds, vox_width, dims);    
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  // we want to find the few branches that are nearest to each voxel
  // possibly we want there to be no maximum distance... which is weird, but more robust I guess.
  // so the best option is to use knn, and match voxel centres to tree segment centres I guess. This has the advantage that
  // it tends not to align leaves to really thick trunks.
  std::vector<int> tree_ids;
  std::vector<int> segment_ids;
  std::vector<std::vector<int> > neighbour_segments; // this looks up into the above two structures
  ForestStructure forest;
 
  { // Tim: this block looks for the closest cylindrical branch segments to each voxel, in order to give the leaves a 'direction' value
    // The reason I use knn (K-nearest neighbour search) is that there is no maximum distance to worry about, and it is fast
    if (!forest.load(trees_file))
    {
      return false;
    }

    size_t num_segments = 0;
    for (auto &tree: forest.trees)
    {
      num_segments += tree.segments().size() - 1;
    }
    size_t num_dense_voxels = 0;
    for (auto &vox: grid.voxels())
    {
      if (vox.density() > 0.0)
        num_dense_voxels++;
    }

    const int search_size = 4; // find the four nearest branch segments. For larger voxels a larger value here would be helpful
    size_t p_size = num_segments;
    size_t q_size = num_dense_voxels;
    Eigen::MatrixXd points_p(3, p_size);
    int i = 0;
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
        points_p.col(i++) = (segment.tip + tree.segments()[segment.parent_id].tip)/2.0;
        tree_ids.push_back(tree_id);
        segment_ids.push_back(segment_id);
      }
    }
    // 2. get 
    neighbour_segments.resize(q_size);
    Eigen::MatrixXd points_q(3, q_size);
    int c = 0;
    for (int i = 0; i<dims[0]; i++)
    {
      for (int j = 0; j<dims[1]; j++)
      {
        for (int k = 0; k<dims[2]; k++)
        {
          int index = grid.getIndex(Eigen::Vector3i(i,j,k));
          double density = grid.voxels()[index].density();
          if (density > 0.0)
          {
            points_q.col(c++) = grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i+0.5, (double)j+0.5, (double)k+0.5);
          }         
        }
      }
    }
    Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    nns->knn(points_q, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
    delete nns;

    // Convert these set of nearest neighbours into surfels
    for (size_t i = 0; i < q_size; i++)
    {
      for (int j = 0; j < search_size && indices(j, i) != Nabo::NNSearchD::InvalidIndex; j++) 
        neighbour_segments[i].push_back(indices(j, i));
    }
  }


  // the density is now stored in grid.voxels()[grid.getIndex(Eigen::Vector3i )].density().
  struct Leaf
  {
    Eigen::Vector3d centre;
    Eigen::Vector3d direction; 
  };
  std::vector<Leaf> leaves;
  int c = 0;
  for (auto &voxel: grid.voxels())
  {
    double leaf_area_per_voxel_volume = voxel.density();
    if (leaf_area_per_voxel_volume <= 0.0)
    {
      continue;
    }
    double desired_leaf_area = leaf_area_per_voxel_volume * vox_width * vox_width * vox_width;
    double num_leaves_d = desired_leaf_area / leaf_area;
    
    // how do we distribute these leaves evenly through the voxel? 
    // we could subvoxel to get a set of point densities...
    // or we could store a scatter matrix per voxel... or just a mean
    // how about we look at the # leaves... if it is more than 8 we randomly place in 2x2x2 subvoxels, if it is more than 27 we randomly place in 3x3x3 subvoxels etc
    int subvoxels = static_cast<int>(std::pow(num_leaves_d, 1.0/3.0));
    double leaves_per_subvoxel = (double)num_leaves_d / (double)(subvoxels * subvoxels * subvoxels);
    double leaves_added = 0;
    double count = 0;
    std::vector<int> &inds = neighbour_segments[c];
    for (int i = 0; i<subvoxels; i++)
    {
      for (int j = 0; j<subvoxels; j++)
      {
        for (int k = 0; k<subvoxels; k++)
        {
          while (leaves_added < leaves_per_subvoxel*count)
          {
            // add leaf
            Eigen::Vector3d minp = grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i, (double)j, (double)k) / (double)subvoxels;
            Eigen::Vector3d maxp = grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i+1, (double)j+1, (double)k+1) / (double)subvoxels;
            Leaf new_leaf;
            new_leaf.centre = Eigen::Vector3d(random(minp[0], maxp[0]), random(minp[1], maxp[1]), random(minp[2], maxp[2]));

            double min_dist = 1e10;
            Eigen::Vector3d closest_point_on_branch(0,0,0);
            int min_ind = 0;
            for (auto &ind: inds)
            {
              auto &tree =  forest.trees[tree_ids[ind]];
              auto &segment = tree.segments()[segment_ids[ind]];
              // get a more accurate distance to each branch segment....
              // e.g. point to branch surface.

              // ADD CODE HERE in order to estimate the closest_point_on_branch
            }

            new_leaf.direction = (new_leaf.centre - closest_point_on_branch).normalized();
            leaves.push_back(new_leaf);
            leaves_added++;
          }
          count++;
        }
      }
    }
    c++;
  }

  // secondly, we place leaves in proportion to the foliage density and (where it is low) we place them to cover real lidar points perhaps
  // better would be to estimate foliage density as a function of direction... in order to distribute the leaf directions
  // perhaps for now we assume that leaves bend downwards, and give them a random range of pitch values...

  // thirdly, we need to save the leaves out.. either as an extension of the tree.txt file, or as an accompanying file. 
  // I wonder if leaf data is enough that a binary format is better?
  // probably we store it as an accompanying data structure for now. It is too early for full integration
  // but we need a function to convert the data to mesh. 
  // maybe we go directly to mesh for the moment.
  Mesh leaf_mesh;
  // could read it from file at this point
  auto &leaf_verts = leaf_mesh.vertices();
  auto &leaf_inds = leaf_mesh.indexList(); // one per triangle, gives the index into the vertices_ array for each corner

  if (leaf_file.empty())
  {
    double droop = 1.0;
    // generate a 2-triangle leaf along y axis
    double len = std::sqrt(leaf_area/2.0);
    leaf_verts.push_back(Eigen::Vector3d(0,0,0));
    leaf_verts.push_back(Eigen::Vector3d(-len/2.0,len,0));
    leaf_verts.push_back(Eigen::Vector3d(len/2.0,len,0));
    leaf_verts.push_back(Eigen::Vector3d(0,2.0*len,-0.5*len * droop));
    leaf_inds.push_back(Eigen::Vector3i(0,1,2));
    leaf_inds.push_back(Eigen::Vector3i(2,1,3));
  }
  else
  {
    readPlyMesh(leaf_file, leaf_mesh);
    // work out its total area:
    double total_area = 0.0;
    for (auto &tri: leaf_inds)
    {
      Eigen::Vector3d side = (leaf_verts[tri[1]] - leaf_verts[tri[0]]).cross(leaf_verts[tri[2]] - leaf_verts[tri[0]]);
      total_area += side.norm() / 2.0;
    }
    double scale = std::sqrt(leaf_area / total_area);
    for (auto &vert: leaf_verts)
    {
      vert *= scale;
    }
  }
      
  
  Mesh mesh;
  auto &verts = mesh.vertices();
  auto &inds = mesh.indexList(); // one per triangle, gives the index into the vertices_ array for each corner   
  
  for (auto &leaf: leaves)
  {
    // 1. convert direction into a transformation matrix...
    Eigen::Matrix3d mat;
    mat.col(0) = leaf.direction;
    mat.col(1) = leaf.direction.cross(Eigen::Vector3d(0,0,1)).normalized();
    mat.col(2) = mat.col(0).cross(mat.col(1));

    int num_verts = (int)verts.size();
    for (auto &tri: leaf_inds)
    {
      inds.push_back(tri + Eigen::Vector3i(num_verts, num_verts, num_verts));
    }
    for (auto &vert: leaf_verts)
    {
      verts.push_back(mat * vert + leaf.centre);
    }
  }          
  writePlyMesh(cloud_stub + "_leaves.ply", mesh);
  return true;
}
}  // namespace ray