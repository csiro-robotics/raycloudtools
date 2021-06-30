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
struct QueueNode
{
  QueueNode(){}
  QueueNode(double distance_to_ground, int index) : distance_to_ground(distance_to_ground), id(index) {}

  double distance_to_ground;
  int id;
};

class QueueNodeComparator 
{ 
public: 
    bool operator() (const QueueNode &p1, const QueueNode &p2) 
    { 
        return p1.distance_to_ground > p2.distance_to_ground; 
    } 
}; 

static const double inf = 1e10;

struct Vertex
{
  Vertex(){}
  Vertex(const Eigen::Vector3d &pos) : pos(pos), parent(-1), branch(-1), distance_to_ground(inf), pilot_id(-1), visited(false) {}
  Eigen::Vector3d pos;
  int parent;
  int branch;
  double distance_to_ground;
  int pilot_id; // if this is a piilot point, then give its id, else it is -1
  bool visited;
};

Trees::Trees(const Cloud &cloud, bool verbose)
{
  std::vector<Vertex> points;  
  points.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
    if (cloud.rayBounded(i))
      points.push_back(Vertex(cloud.ends[i]));
  if (verbose)
  {
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
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
    points_p.col(i) = points[i].pos;
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
  for (unsigned int i = 0; i<points.size(); i++)
  {
    int num_neighbours;
    for (num_neighbours = 0; num_neighbours < search_size && indices(num_neighbours, i) > -1; num_neighbours++);
    Eigen::Vector3d &point = points[i].pos;
    bool any_below = false;
    for (int j = 0; j<num_neighbours && !any_below; j++)
    {
      int id = indices(j, i);
      if (points[id].pos[2] < point[2])
        any_below = true;
    }
    if (!any_below)
    {
      if (points[i].pos[2] < 0.05)
        ground_points.push_back(i);
    }
  }
  
  // 2b. climb up from lowest points...
	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;

  for (auto &ground_id: ground_points)
  {
    points[ground_id].distance_to_ground = 0;
    closest_node.push(QueueNode(0, ground_id));
  }

	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!points[node.id].visited)
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        if (node.distance_to_ground + dists(i, node.id) < points[child].distance_to_ground)
        {
					points[child].distance_to_ground = node.distance_to_ground + dists(i, node.id);
          points[child].parent = node.id;
          points[child].branch = points[node.id].branch;
          points[child].pilot_id = points[node.id].pilot_id;
					closest_node.push(QueueNode(points[child].distance_to_ground, child));
				}
			}
		  points[node.id].visited = true;
		}
	}
/*
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

  TreesNode root;
  root.roots = ground_points;
  tree_nodes.push_back(root);

  // now trace from root tree nodes upwards, getting node centroids
  for (unsigned int tn = 0; tn < tree_nodes.size(); tn++)
  {
    for (int i = 0; i<vertices; i++)
      visited[i] = false;   
    double thickness = node_separation;
    if (tree_nodes[tn].parent >= 0)
      thickness = 2.0*tree_nodes[tree_nodes[tn].parent].radius;
    double max_dist_from_ground = tree_nodes[tn].min_dist_from_ground + thickness;
    Eigen::Vector3d centroid_sqr(0,0,0);
    std::vector<int> child_roots;
    std::vector<int> nodes = tree_nodes[tn].roots;
    for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
    {
      int i = nodes[ijk];
      if (visited[i])
        continue;
      visited[i] = true;
      tree_nodes[tn].centroid += points[i];
      centroid_sqr += Eigen::Vector3d(sqr(points[i][0]), sqr(points[i][1]), sqr(points[i][2]));
      tree_nodes[tn].num_points++;
      int id = i;
      for (auto &child: children[id])
      {
        if (visited[child])
          continue;
        if (distances_to_ground[child] < max_dist_from_ground) // in same slot, so accumulate
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
      centroid_sqr /= (double)tree_nodes[tn].num_points;
    }
    Eigen::Vector3d centroid2(sqr(tree_nodes[tn].centroid[0]), sqr(tree_nodes[tn].centroid[1]), sqr(tree_nodes[tn].centroid[2]));
    Eigen::Vector3d variance = centroid_sqr - centroid2;
    tree_nodes[tn].radius = std::pow(variance[0] * variance[1] * variance[2], 1.0/6.0);
    if (tree_nodes[tn].radius > 1.2)
      tree_nodes[tn].radius = node_separation;
    tree_nodes[tn].radius = std::max(tree_nodes[tn].radius, 0.0125);
  //  std::cout << "cent: " << centroid_sqr.transpose() << ", cen2: " << centroid2.transpose() << " = rad: " << tree_nodes[tn].radius << std::endl;
    
    // floodfill graph to find separate cliques
    for (int i = 0; i<vertices; i++)
      visited[i] = false;   
    int num_added = 0;
    for (auto &root: child_roots)
    {
      if (visited[root])
        continue;
      TreesNode new_node;
      new_node.parent = tn;
      new_node.roots.push_back(root);
      new_node.min_dist_from_ground = max_dist_from_ground;
      for (size_t ijk = 0; ijk<new_node.roots.size(); ijk++)
      {
        int i = new_node.roots[ijk];
        visited[i] = true;
        for (auto &j: child_roots)
        {
          if (visited[j])
            continue;
          if ((points[i] - points[j]).norm() < tree_nodes[tn].radius)
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
        starts.push_back(tree_nodes[tree_node.parent].centroid);
        ends.push_back(tree_node.centroid);
  //      radii.push_back(tree_node.radius);
      }
    }
    DebugDraw::instance()->drawLines(starts, ends);
//    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
  }

  // generate children links
  for (unsigned int i = 0; i<tree_nodes.size(); i++)
  {
    if (tree_nodes[i].parent >= 0)
    {
      tree_nodes[tree_nodes[i].parent].children.push_back(i);
    }
    else if (tree_nodes[i].parent == -2)
      root_nodes.push_back(i);
  }
  // generate local parent links
  for (auto &root: root_nodes)
  {
    int child_id = 0;
    tree_nodes[root].id = child_id++;
    std::vector<int> children = tree_nodes[root].children;
    for (unsigned int c = 0; c<children.size(); c++)
    {
      tree_nodes[children[c]].id = child_id++;
      for (auto i: tree_nodes[children[c]].children)
        children.push_back(i);
    }
  }*/
}

bool Trees::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree structure. One row per tree, which repeats: 'x,y,z,radius,parent_id, ' per segment" << std::endl;
  for (auto &root: root_nodes)
  {
    ofs << tree_nodes[root].centroid[0] << "," << tree_nodes[root].centroid[1] << "," << tree_nodes[root].centroid[2] << "," << tree_nodes[root].radius << ",-1";

    std::vector<int> children = tree_nodes[root].children;
    for (unsigned int c = 0; c<children.size(); c++)
    {
      TreesNode &node = tree_nodes[children[c]];
      ofs << ", " << node.centroid[0] << "," << node.centroid[1] << "," << node.centroid[2] << "," << node.radius << "," << tree_nodes[node.parent].id;
      for (auto i: tree_nodes[children[c]].children)
        children.push_back(i);
    }
    ofs << std::endl;
  }
  return true;
}

} // namespace ray
