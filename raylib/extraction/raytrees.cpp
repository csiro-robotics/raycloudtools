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
  Vertex(const Eigen::Vector3d &pos) : pos(pos), parent(-1), distance_to_ground(inf), radius(0.0), visited(false) {}
  Eigen::Vector3d pos;
  int parent;
  double distance_to_ground;
  double radius;
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
 //   DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
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
  double start_radius = 0.16;
  for (auto &ground_id: ground_points)
  {
    points[ground_id].distance_to_ground = 0;
    points[ground_id].radius = start_radius;
    closest_node.push(QueueNode(0, ground_id));
  }

	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!points[node.id].visited)
    {
      // just estimate radius here, based on parent's radius estimate
      double min_dist = 2.0*points[node.id].radius;
      double min_dist_sqr = sqr(min_dist);
      Eigen::Vector3d sum = points[node.id].pos;
      Eigen::Matrix3d sum_sqr = points[node.id].pos * points[node.id].pos.transpose();
      double num = 1;
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int neighbour = indices(i, node.id);
        if (points[neighbour].distance_to_ground > points[node.id].distance_to_ground - min_dist)
        {
          if ((points[node.id].pos - points[neighbour].pos).squaredNorm() < min_dist_sqr)
          {
            sum += points[neighbour].pos;
            sum_sqr += points[neighbour].pos * points[neighbour].pos.transpose();
            num++;
          }
        }
      }
      sum /= num;
      Eigen::Matrix3d scatter = (sum_sqr/num) - sum*sum.transpose(); // TODO: fix, this can't be right
      double radius = pow(std::abs(scatter.determinant()), 1.0/6.0);

      if (num > 3 && points[node.id].radius > radius)
      {
        points[node.id].radius += (radius - points[node.id].radius) * 0.5;
      }

      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        if (node.distance_to_ground + dists(i, node.id) < points[child].distance_to_ground)
        {
					points[child].distance_to_ground = node.distance_to_ground + dists(i, node.id);
          points[child].parent = node.id;
          points[child].radius = points[node.id].radius;
					closest_node.push(QueueNode(points[child].distance_to_ground, child));
				}
			}
		  points[node.id].visited = true;
		}
	}
  if (verbose)
  {
    std::vector<Eigen::Vector3d> ps(points.size());
    std::vector<double> vs(points.size());
    for (size_t i = 0; i<points.size(); i++)
    {
      ps[i] = points[i].pos;
      vs[i] = points[i].radius / start_radius;
    }
    DebugDraw::instance()->drawCloud(ps, vs, 0);
  }  

  const double node_separation = 0.16;
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<Eigen::Vector3d> colours;
    for (size_t i = 0; i<points.size(); i++)
    {
      if (points[i].parent != -1)
      {
        starts.push_back(points[points[i].parent].pos);
        ends.push_back(points[i].pos);
        int slot = (int)(points[i].distance_to_ground / node_separation);
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
  std::vector< std::vector<int> > children(points.size());
  for (size_t i = 0; i<points.size(); i++)
  {
    if (points[i].parent != -1)
      children[points[i].parent].push_back((int)i);
  }
  // now generate initial sections
  for (size_t i = 0; i<points.size(); i++)
    points[i].visited = false;

  TreesNode root;
  root.roots = ground_points;
  tree_nodes.push_back(root);

  // now trace from root tree nodes upwards, getting node centroids
  for (size_t tn = 0; tn < tree_nodes.size(); tn++)
  {
    for (size_t i = 0; i<points.size(); i++)
      points[i].visited = false;   
    double thickness = node_separation;
    if (tree_nodes[tn].parent >= 0)
      thickness = 4.0*tree_nodes[tree_nodes[tn].parent].radius;
    double max_dist_from_ground = tree_nodes[tn].min_dist_from_ground + thickness;
    Eigen::Vector3d centroid_sqr(0,0,0);
    std::vector<int> child_roots;
    std::vector<int> nodes = tree_nodes[tn].roots;
    std::vector<Eigen::Vector3d> radius_points;
    for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
    {
      int i = nodes[ijk];
      if (points[i].visited)
        continue;
      points[i].visited = true;
      tree_nodes[tn].centroid += points[i].pos;
      tree_nodes[tn].num_points++;
      radius_points.push_back(points[i].pos);
      int id = i;
      for (auto &child: children[id])
      {
        if (points[child].visited)
          continue;
        if (points[child].distance_to_ground < max_dist_from_ground) // in same slot, so accumulate
        {
          nodes.push_back(child); // so we recurse on this child too
        }
        else // not same slot, so add to child roots 
        {
          child_roots.push_back(child);
          points[child].visited = true;
        }
      }
    }
    if (tree_nodes[tn].num_points != (int)nodes.size())
      std::cout << "node size isn't a proxy for num points" << std::endl;
    if (tree_nodes[tn].num_points > 0)
      tree_nodes[tn].centroid /= (double)tree_nodes[tn].num_points;
    
    Eigen::Vector3d dir(0,0,1);
    if (tree_nodes[tn].parent != -1)
      dir = (tree_nodes[tn].centroid - tree_nodes[tree_nodes[tn].parent].centroid).normalized();
    double rad_sqr = 0.0;
    for (auto &node: nodes)
    {
      Eigen::Vector3d offset = points[node].pos - tree_nodes[tn].centroid;
      Eigen::Vector3d p = offset - dir*offset.dot(dir);
      rad_sqr += p.squaredNorm();
    }
    if (tree_nodes[tn].num_points > 4 || tree_nodes[tn].parent == -1)
    {
      rad_sqr /= (double)tree_nodes[tn].num_points;
      tree_nodes[tn].radius = std::sqrt(rad_sqr);
    }
    else
    {
      tree_nodes[tn].radius = 1.2*tree_nodes[tree_nodes[tn].parent].radius;
    }

    {
 //     DebugDraw::instance()->drawCloud(radius_points, 0.5, 0); 
 //     std::cout << "dir: " << dir.transpose() << " radius: " << tree_nodes[tn].radius << " of " << tree_nodes[tn].num_points << " points" << std::endl;
    }

    if (tree_nodes[tn].radius > 1.0)
      tree_nodes[tn].radius = 0.2;//node_separation;



    tree_nodes[tn].radius = std::max(tree_nodes[tn].radius, 0.0125);
    if (tree_nodes[tn].parent != -1)
      tree_nodes[tn].radius = std::min(tree_nodes[tn].radius, 1.1*tree_nodes[tree_nodes[tn].parent].radius);
      
  //  std::cout << "cent: " << centroid_sqr.transpose() << ", cen2: " << centroid2.transpose() << " = rad: " << tree_nodes[tn].radius << std::endl;
    
    // floodfill graph to find separate cliques
    for (size_t i = 0; i<points.size(); i++)
      points[i].visited = false;   
    for (auto &root: child_roots)
    {
      if (points[root].visited)
        continue;
      TreesNode new_node;
      new_node.parent = (int)tn;
      new_node.roots.push_back(root);
      new_node.min_dist_from_ground = max_dist_from_ground;
      for (size_t ijk = 0; ijk<new_node.roots.size(); ijk++)
      {
        int i = new_node.roots[ijk];
        points[i].visited = true;
        for (auto &j: child_roots)
        {
          if (points[j].visited)
            continue;
          if ((points[i].pos - points[j].pos).norm() < 1.5*tree_nodes[tn].radius)
          {
            new_node.roots.push_back(j);
            points[j].visited = true;
          }
        }        
      }
      tree_nodes.push_back(new_node);
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
    std::vector<double> radii;
    for (auto &tree_node: tree_nodes)
    {
      if (tree_node.parent >= 0)
      {
        if ((tree_nodes[tree_node.parent].centroid - tree_node.centroid).norm() < 0.01)
          continue;
        starts.push_back(tree_nodes[tree_node.parent].centroid);
        ends.push_back(tree_node.centroid);
        radii.push_back(tree_node.radius);
      }
    }
    DebugDraw::instance()->drawLines(starts, ends);
    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
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
