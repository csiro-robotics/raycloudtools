// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raysegment.h"
#include "rayterrain.h"
#include <nabo/nabo.h>

namespace ray
{

void connectPointsShortestPath(std::vector<Vertex> &points, std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node, double distance_limit, double gravity_factor)
{
  // 1. get nearest neighbours
  const int search_size = 20;
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++) 
    points_p.col(i) = points[i].pos;
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, distance_limit);
  
  // 2b. climb up from lowest points...
	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!points[node.id].visited)
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        double dist2 = dists2(i, node.id);
        double dist = std::sqrt(dist2);
        double new_score = 0;
        Eigen::Vector3d dif = (points[child].pos - points[node.id].pos).normalized();
        Eigen::Vector3d dir(0,0,1);
        int ppar = points[node.id].parent;
        if (ppar != -1)
        {
          if (points[ppar].parent != -1)
            dir = (points[node.id].pos - points[points[ppar].parent].pos).normalized(); // this is a bit smoother than...
          else
            dir = (points[node.id].pos - points[ppar].pos).normalized();  // ..just this
        }
        const double power = 2.0;
        dist2 /= std::pow(std::max(0.001, dif.dot(dir)), power);

        if (gravity_factor > 0.0)
        {
          Eigen::Vector3d to_node = points[node.id].pos - points[node.root].pos;
          to_node[2] = 0.0;
          double lateral_sqr = to_node.squaredNorm();
          double gravity_scale = 1.0 + gravity_factor*lateral_sqr; // the squaring means gravity plays little role for normal trees, kicking in stronger on outlier lateral ones
          dist2 *= gravity_scale;
        }

        dist2 /= node.radius;
        new_score = node.score + dist2;
        if (new_score < points[child].score)
        {
					points[child].score = new_score;
					points[child].distance_to_ground = node.distance_to_ground + dist;
          points[child].parent = node.id;
          points[child].root = node.root;
					closest_node.push(QueueNode(points[child].distance_to_ground, points[child].score, node.radius, node.root, child));
				}
			}
		  points[node.id].visited = true;
		}
	}
}

std::vector< std::vector<int> > getRootsAndSegment(std::vector<Vertex> &points, Cloud &cloud, const Mesh &mesh, double max_diameter, double distance_limit, double height_min, double gravity_factor)
{
  points.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    if (cloud.rayBounded(i))
    {
      points.push_back(Vertex(cloud.ends[i]));
    }      
  }

  const double pixel_width = max_diameter;
  Eigen::Vector3d box_min, box_max;
  cloud.calcBounds(&box_min, &box_max);
	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;

  int roots_start = (int)points.size();
  for (auto &vert: mesh.vertices())
  {
    points.push_back(Vertex(vert));
  }
  Eigen::ArrayXXd lowfield;
  mesh.toHeightField(lowfield, box_min, box_max, pixel_width);

  // set the height field to the height of the canopy above the ground
  Eigen::ArrayXXd heightfield = Eigen::ArrayXXd::Constant((int)lowfield.rows(), (int)lowfield.cols(), -1e10);
  for (auto &point: points)
  {
    Eigen::Vector3i index = ((point.pos - box_min)/pixel_width).cast<int>();
    heightfield(index[0], index[1]) = std::max(heightfield(index[0], index[1]), point.pos[2]);
  }
  for (int i = 0; i<heightfield.rows(); i++)
  {
    for (int j = 0; j<heightfield.cols(); j++)
    {
      heightfield(i,j) = std::max(1e-10, heightfield(i,j) - lowfield(i,j));
    }
  }

  for (int ind = roots_start; ind<(int)points.size(); ind++) 
  {
    points[ind].distance_to_ground = 0.0;
    points[ind].score = 0.0;
    points[ind].root = ind;
    Eigen::Vector3i index = ((points[ind].pos - box_min)/pixel_width).cast<int>();
    closest_node.push(QueueNode(0, 0, heightfield(index[0], index[1]), ind, ind));
  }

  connectPointsShortestPath(points, closest_node, distance_limit, gravity_factor);

  // now create a count array
  Eigen::ArrayXXi counts = Eigen::ArrayXXi::Constant((int)heightfield.rows(), (int)heightfield.cols(), 0);
  Eigen::ArrayXXd heights = Eigen::ArrayXXd::Constant((int)heightfield.rows(), (int)heightfield.cols(), 0);
  for (auto &point: points)
  {
    if (point.root == -1)
      continue;
    Eigen::Vector3i index = ((points[point.root].pos - box_min)/pixel_width).cast<int>();    
    counts(index[0], index[1])++;
    heights(index[0], index[1]) = std::max(heights(index[0], index[1]), point.pos[2] - lowfield(index[0], index[1]));
  }
  // now create a 2x2 summed array:
  Eigen::ArrayXXi sums = Eigen::ArrayXXi::Constant((int)counts.rows(), (int)counts.cols(), 0);
  for (int i = 0; i<(int)sums.rows(); i++)
  {
    for (int j = 0; j<(int)sums.cols(); j++)
    {
      int i2 = std::min(i+1, (int)sums.rows()-1);
      int j2 = std::min(j+1, (int)sums.cols()-1);
      sums(i,j) = counts(i,j)+counts(i,j2)+counts(i2,j)+counts(i2,j2);
    }
  }
  // now find the best 2x2 sum for each cell:
  std::vector<Eigen::Vector2i> bests((int)counts.rows() * (int)counts.cols());
  for (int x = 0; x<(int)sums.rows(); x++)
  {
    for (int y = 0; y<(int)sums.cols(); y++)
    {
      Eigen::Vector2i best_index(-1,-1);
      int largest_sum = -1;
      for (int i = std::max(0, x-1); i<=x; i++)
      {
        for (int j = std::max(0, y-1); j<=y; j++)
        {
          if (sums(i,j) > largest_sum)
          {
            largest_sum = sums(i,j);
            best_index = Eigen::Vector2i(i,j);
          }
        }
      }
      bests[x + sums.rows()*y] = best_index;
    }
  }
  // next we need to find the highest point for each cell....
  Eigen::ArrayXXd max_heights = Eigen::ArrayXXd::Constant((int)counts.rows(), (int)counts.cols(), 0);
  for (int i = 0; i<(int)sums.rows(); i++)
  {
    for (int j = 0; j<(int)sums.cols(); j++)
    {
      Eigen::Vector2i best_index = bests[i + sums.rows()*j];
      double max_height = 0.0;    
      for (int x = best_index[0]; x<std::min(best_index[0]+2, (int)sums.rows()); x++)
      {
        for (int y = best_index[1]; y<std::min(best_index[1]+2, (int)sums.cols()); y++)
        {
          if (bests[x + (int)sums.rows()*y] == best_index)
          {
            max_height = std::max(max_height, heights(x,y));
          }
        }
      }
      max_heights(best_index[0], best_index[1]) = max_height;
    }
  }  

  std::vector< std::vector<int> > roots_lists(sums.rows() * sums.cols());
  for (int i = roots_start; i < (int)points.size(); i++)
  {
    Eigen::Vector3i index = ((points[i].pos - box_min)/pixel_width).cast<int>();    
    Eigen::Vector2i best_index = bests[index[0] + (int)sums.rows()*index[1]];
    double max_height = max_heights(best_index[0], best_index[1]);
    if (max_height >= height_min)
    {
      int id = best_index[0] + (int)sums.rows()*best_index[1];
      roots_lists[id].push_back(i);
    }
  }
  std::vector< std::vector<int> > roots_set; // contiguous form of root_lists
  for (int i = 0; i<(int)sums.rows(); i++)
  {
    for (int j = 0; j<(int)sums.cols(); j++)
    {
      auto &roots = roots_lists[i + (int)sums.rows()*j];
      if (roots.size() > 0)
      {
        roots_set.push_back(roots);
      }
    }
  }

  return roots_set;
}

} // namespace ray
