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
void connectPointsShortestPath(std::vector<Vertex> &points, std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node, double distance_limit)
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
        double dist = std::sqrt(dists2(i, node.id));
        double new_dist = node.distance_to_ground + dist/node.radius;
        double new_score = 0;
        #if defined MINIMISE_SQUARE_DISTANCE
        dist *= dist;
        #endif
        #if defined MINIMISE_ANGLE
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
        dist /= std::pow(std::max(0.001, dif.dot(dir)), power);
        #endif
        dist /= node.radius;
        #if defined MINIMISE_SQUARE_DISTANCE || defined MINIMISE_ANGLE
        new_score = node.score + dist;
        if (new_score < points[child].score)
        #else
        if (new_dist < points[child].distance_to_ground)
        #endif
        {
					points[child].score = new_score;
					points[child].distance_to_ground = new_dist;
          points[child].parent = node.id;
          points[child].root = node.root;
					closest_node.push(QueueNode(points[child].distance_to_ground, points[child].score, node.radius, node.root, child));
				}
			}
		  points[node.id].visited = true;
		}
	}
}

void segmentTrees(Cloud &cloud, double max_diameter, double gradient, double distance_limit, double height_min)
{
  std::cout << "segmenting cloud with tree diameter: " << max_diameter << " and gradient: " << gradient << std::endl;
  std::vector<Vertex> points;  
  points.reserve(cloud.ends.size());
  std::vector<Eigen::Vector3d> raw_points;
  raw_points.reserve(cloud.ends.size());
  std::vector<int> original_ids;
  original_ids.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    if (cloud.rayBounded(i))
    {
      points.push_back(Vertex(cloud.ends[i]));
      raw_points.push_back(cloud.ends[i]);
      original_ids.push_back(i);
    }      
  }

  const double pixel_width = max_diameter;
  Eigen::Vector3d box_min, box_max;
  cloud.calcBounds(&box_min, &box_max);
	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;

  Terrain hull;
  double pix_width = 2.0 * cloud.estimatePointSpacing();
  hull.growUpwardsFast(raw_points, pix_width, box_min, box_max, gradient);
  int roots_start = (int)points.size();
  for (auto &vert: hull.mesh().vertices())
  {
    points.push_back(Vertex(vert));
  }
  Eigen::ArrayXXd lowfield;
  hull.mesh().toHeightField(lowfield, box_min, box_max, pixel_width);

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

  connectPointsShortestPath(points, closest_node, distance_limit);

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
          if (bests[x + sums.rows()*y] == best_index)
          {
            max_height = std::max(max_height, heights(x,y));
          }
        }
      }
      max_heights(i, j) = max_height;
    }
  }  

  // now colour the cloud based on start index
  for (int i = 0; i<roots_start; i++)
  {
    RGBA &colour = cloud.colours[original_ids[i]];
    Vertex &point = points[i];
    if (point.root == -1)
    {
      colour.red = colour.green = colour.blue = 0;
      continue;
    }
    Eigen::Vector3i index = ((points[point.root].pos - box_min)/pixel_width).cast<int>();    
    Eigen::Vector2i best_index = bests[index[0] + sums.rows()*index[1]];
    double max_height = max_heights(index[0], index[1]);
    if (max_height < height_min)
    {
      colour.red = colour.green = colour.blue = 0;
      continue;
    }
    index[0] = best_index[0]; index[1] = best_index[1];
    srand(index[0] + 3127*index[1]);
    colour.red   = uint8_t(50 + rand()%205);
    colour.green = uint8_t(50 + rand()%205);
    colour.blue  = uint8_t(50 + rand()%205);
  }
}

} // namespace ray
