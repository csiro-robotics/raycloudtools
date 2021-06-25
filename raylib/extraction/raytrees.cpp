// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include "../raydebugdraw.h"

#include <nabo/nabo.h>
#include <queue>

namespace ray
{
struct Node
{
  Node(){}
  Node(double distance_to_ground, int index) : distance_to_ground(distance_to_ground), id(index) {}

  double distance_to_ground;
  int id;
};

class myComparator 
{ 
public: 
    bool operator() (const Node &p1, const Node &p2) 
    { 
        return p1.distance_to_ground > p2.distance_to_ground; 
    } 
}; 

Trees::Trees(const Cloud &cloud, bool verbose)
{
  std::vector<Eigen::Vector3d> points;  
  points.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
    if (cloud.rayBounded(i))
      points.push_back(cloud.ends[i]);
  if (verbose)
  {
    DebugDraw::instance()->drawCloud(points, 0.5, 0);
  }

  // OK, the planned algorithm is as follows:
  // 1. get nearest neighbours for all points
  // 2. walk from the lowest points upwards, for each point above: find its shortest path to a root point, to develop a set of trees
  // 3. floodfill neighbours according to quantised distance from ground
  // 4. generate skeletons from connecting these groups together
  // 5. calculate radius of skeleton at each node
  // 6. try to remove clutter
  // 7. render the trees

  // 1. get nearest neighbours
  const int search_size = 20;
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++) 
    points_p.col(i) = points[i];
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
  Eigen::MatrixXd dists = dists2;
  for (int i = 0; i<dists2.rows(); i++)
  {
    for (int j = 0; j<dists2.cols(); j++)
      dists(i, j) = std::sqrt(dists(i,j));
  }

  // 2. walk the lowest points upwards....
  // 2a. get the lowest points
  std::vector<int> ground_points;
 // std::vector<int> visited(points.size());
 // for (auto &b: visited)
 //   b = 0;
  for (unsigned int i = 0; i<points.size(); i++)
  {
    int num_neighbours;
    for (num_neighbours = 0; num_neighbours < search_size && indices(num_neighbours, i) > -1; num_neighbours++);
    Eigen::Vector3d &point = points[i];
    bool any_below = false;
    for (int j = 0; j<num_neighbours && !any_below; j++)
    {
      int id = indices(j, i);
      if (points[id][2] < point[2])
        any_below = true;
    }
    if (!any_below)
    {
      if (points[i][2] < 0.05)
        ground_points.push_back(i);
 //     if (points[i][2] > 0.05)
 //       std::cout << "spurious: " << points[i].transpose() << std::endl;
  //    visited[i] = 1;
    }
  }
  
  // 2b. climb up from lowest points...
  const int vertices = (int)points.size();
  const double inf = 1e10;

	std::priority_queue<Node, std::vector<Node>, myComparator> closest_node;

	std::vector<bool> visited(vertices);
	std::vector<double> distances_to_ground(vertices);
	std::vector<int> parent_id(vertices);
	for (int i = 0; i < vertices; i++)
	{
		visited[i] = 0;
		distances_to_ground[i] = inf;
    parent_id[i] = -1;
	}

  for (auto &ground_id: ground_points)
 // for (int i = 0; i<1; i++) // int ground_id = 0;
  {
  //  int ground_id = ground_points[i];
    distances_to_ground[ground_id] = 0;
    parent_id[ground_id] = -1;
    closest_node.push(Node(0, ground_id));
  }

	while(!closest_node.empty())
  {
		Node node = closest_node.top(); closest_node.pop();
		if(!visited[node.id])
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        if (node.distance_to_ground + dists(i, node.id) < distances_to_ground[child])
        {
					distances_to_ground[child] = node.distance_to_ground + dists(i, node.id);
          parent_id[child] = node.id;
					closest_node.push(Node(distances_to_ground[child], child));
				}
			}
		  visited[node.id] = true;
		}
	}
  const double node_separation = 0.16;
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<Eigen::Vector3d> colours;
    for (int i = 0; i<vertices; i++)
    {
      if (parent_id[i] != -1)
      {
        starts.push_back(points[parent_id[i]]);
        ends.push_back(points[i]);
        int slot = (int)(distances_to_ground[i] / node_separation);
        srand(slot);
        Eigen::Vector3d col;
        col[0] = (double)(rand()%1000)/1000.0;
        col[1] = (double)(rand()%1000)/1000.0;
        col[2] = (double)(rand()%1000)/1000.0;
        colours.push_back(col);
      }
    }
    DebugDraw::instance()->drawLines(starts, ends, colours);
  }

  // 2c. generate skeletons. 
  // First generate children
  std::vector< std::vector<int> > children(vertices);
  for (int i = 0; i<vertices; i++)
  {
    if (parent_id[i] != -1)
      children[parent_id[i]].push_back(i);
  }
  // now generate initial sections
  for (int i = 0; i<vertices; i++)
    visited[i] = false;
  struct TreeNode
  {
    TreeNode() : centroid(0,0,0), centroid_sqr(0,0,0), num_points(0), parent(-1) {}
    Eigen::Vector3d centroid;
    Eigen::Vector3d centroid_sqr;
    int num_points;
    int parent;
    std::vector<int> roots; // root points
  };    
  std::vector<TreeNode> tree_nodes;

  // how do we define separate paths.......
  // flood fill the ground points by slot (quantisation of distance from ground)
 /* for (auto &root: ground_points)
  {
    if (visited[root])
      continue;
    int slot = (int)(distances_to_ground[root] / node_separation);
    if (slot > 0)
      continue;
    TreeNode new_node;
    new_node.roots.push_back(root);
    for (size_t ijk = 0; ijk<new_node.roots.size(); ijk++)
    {
      int i = new_node.roots[ijk];
      visited[i] = true;
      for (int j = 0; j<search_size && indices(j, i) > -1; j++)
      {
        int neighbour_id = indices(j, i);
        if (visited[neighbour_id])
          continue;
        int neighbour_slot = (int)(distances_to_ground[neighbour_id] / node_separation);
        if (neighbour_slot == slot)
        {
          new_node.roots.push_back(neighbour_id);
          visited[neighbour_id] = true;
        }
      }    
    }
    tree_nodes.push_back(new_node);
  }*/

  TreeNode root;
  root.roots = ground_points;
  tree_nodes.push_back(root);

  // now trace from root tree nodes upwards, getting node centroids
  for (unsigned int tn = 0; tn < tree_nodes.size(); tn++)
  {
    for (int i = 0; i<vertices; i++)
      visited[i] = false;   
    std::vector<int> child_roots;
    std::vector<int> nodes = tree_nodes[tn].roots;
    for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
    {
      int i = nodes[ijk];
      if (visited[i])
        continue;
      int slot = (int)(distances_to_ground[i]/node_separation);
      visited[i] = true;
      tree_nodes[tn].centroid += points[i];
      tree_nodes[tn].centroid_sqr += Eigen::Vector3d(sqr(points[i][0]), sqr(points[i][1]), sqr(points[i][2]));
      tree_nodes[tn].num_points++;
      int id = i;
      for (auto &child: children[id])
      {
        if (visited[child])
          continue;
        int child_slot = (int)(distances_to_ground[child]/node_separation);
        if (child_slot == slot) // in same slot, so accumulate
        {
          nodes.push_back(child); // so we recurse on this child too
        }
        else // not same slot, so add to child roots 
        {
          child_roots.push_back(child);
          visited[child] = true;
        }
      }
    }
    if (tree_nodes[tn].num_points > 0)
    {
      tree_nodes[tn].centroid /= (double)tree_nodes[tn].num_points;
      tree_nodes[tn].centroid_sqr /= (double)tree_nodes[tn].num_points;
    }
    
    // floodfill graph to find separate cliques
    for (int i = 0; i<vertices; i++)
      visited[i] = false;   
    int num_added = 0;
    for (auto &root: child_roots)
    {
      if (visited[root])
        continue;
      TreeNode new_node;
      new_node.parent = tn;
      new_node.roots.push_back(root);
      for (size_t ijk = 0; ijk<new_node.roots.size(); ijk++)
      {
        int i = new_node.roots[ijk];
        visited[i] = true;
        for (auto &j: child_roots)
        {
          if (visited[j])
            continue;
          if ((points[i] - points[j]).norm() < node_separation)
          {
            new_node.roots.push_back(j);
            visited[j] = true;
          }
        }        
      }
      tree_nodes.push_back(new_node);
      num_added++;
    }
  }
  // remove node 0
  for (auto &tree_node: tree_nodes)
  {
    if (tree_node.parent >= 0 && tree_nodes[tree_node.parent].parent == -1)
      tree_node.parent = -2;
  }

  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
 //   std::vector<double> radii;
    for (auto &tree_node: tree_nodes)
    {
      if (tree_node.parent >= 0)
      {
        if ((tree_nodes[tree_node.parent].centroid - tree_node.centroid).norm() < 0.01)
          continue;
        Eigen::Vector3d sigma = tree_node.centroid_sqr - Eigen::Vector3d(sqr(tree_node.centroid[0]), sqr(tree_node.centroid[1]), sqr(tree_node.centroid[2]));
        if (sigma.norm() < 0.01)
          continue;
        starts.push_back(tree_nodes[tree_node.parent].centroid);
        ends.push_back(tree_node.centroid);
/*        sigma[0] = std::sqrt(sigma[0]);
        sigma[1] = std::sqrt(sigma[1]);
        sigma[2] = std::sqrt(sigma[2]);
        radii.push_back(sigma.norm());*/
      }
    }
    DebugDraw::instance()->drawLines(starts, ends);
//    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
  }
}


bool Trees::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree base location list: x, y, z, radius" << std::endl;
/*  for (auto &trunk: trunk_bases)
  {
    Eigen::Vector3d base = trunk.centre - vector3d(trunk.lean, 1)*trunk.length*0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << trunk.radius << std::endl;
  }*/
  return true;
}

} // namespace ray
