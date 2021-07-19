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
  QueueNode(double distance_to_ground, double score, int index) : distance_to_ground(distance_to_ground), score(score), id(index) {}

  double distance_to_ground;
  double score;
  int id;
};

//#define MINIMISE_SCORE // currently the square of distance
class QueueNodeComparator 
{ 
public: 
    bool operator() (const QueueNode &p1, const QueueNode &p2) 
    { 
#if defined MINIMISE_SCORE
        return p1.score > p2.score; 
#else
        return p1.distance_to_ground > p2.distance_to_ground; 
#endif
    } 
}; 

static const double inf = 1e10;

struct Vertex
{
  Vertex(){}
  Vertex(const Eigen::Vector3d &pos) : pos(pos), parent(-1), distance_to_ground(inf), score(inf), visited(false) {}
  Eigen::Vector3d pos;
  int parent;
  double distance_to_ground;
  double score;
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
      if (points[i].pos[2] < 0.15)
        ground_points.push_back(i);
    }
  }

  std::cout << "num ground points: " << ground_points.size() << std::endl;
  
  // 2b. climb up from lowest points...
	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;
  for (auto &ground_id: ground_points)
  {
    points[ground_id].distance_to_ground = 0;
    points[ground_id].score = 0;
    closest_node.push(QueueNode(0, 0, ground_id));
  }

	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!points[node.id].visited)
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        #if defined MINIMISE_SCORE
        if (node.score + dists2(i, node.id) < points[child].score)
        #else
        if (node.distance_to_ground + dists(i, node.id) < points[child].distance_to_ground)
        #endif
        {
					points[child].score = node.score + dists2(i, node.id);
					points[child].distance_to_ground = node.distance_to_ground + dists(i, node.id);
          points[child].parent = node.id;
					closest_node.push(QueueNode(points[child].distance_to_ground, points[child].score, child));
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
      vs[i] = 1.0;
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
  root.radius = 0.1;
  tree_nodes.push_back(root);

  // now trace from root tree nodes upwards, getting node centroids
  // a tree_node is a segment, and these are added as we iterate through the list
  for (size_t tn = 0; tn < tree_nodes.size(); tn++)
  {
    for (size_t i = 0; i<points.size(); i++)
      points[i].visited = false;   
    double thickness = node_separation;
    if (tree_nodes[tn].parent >= 0)
      thickness = 4.0*tree_nodes[tree_nodes[tn].parent].radius;
    double max_dist_from_ground = tree_nodes[tn].min_dist_from_ground + thickness;

    std::vector<int> nodes = tree_nodes[tn].roots;
    std::cout << "tree " << tn << " roots: " << tree_nodes[tn].roots.size() << ", ends: " << tree_nodes[tn].ends.size() << std::endl;
    bool extract_from_ends = tree_nodes[tn].ends.size() > 0;
    if (!extract_from_ends)
    {
      // 1. find all the points in this tree node:
      for (size_t i = 0; i<points.size(); i++)
        points[i].visited = false;  
      for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
      {
        int i = nodes[ijk];
        if (points[i].visited)
          continue;
        points[i].visited = true;
        for (auto &child: children[i])
        {
          if (points[child].visited)
            continue;
          if (points[child].distance_to_ground < max_dist_from_ground) // in same slot, so accumulate
            nodes.push_back(child); // so we recurse on this child too
          else 
          {
            tree_nodes[tn].ends.push_back(child); 
            points[child].visited = true;
          }
        }
      }
      if (verbose)
        std::cout << "no ends, so found " << tree_nodes[tn].ends.size() << " ends" << std::endl;
    }
    // TODO: what if we have found 0 ends here? i.e. the end of the branch...?
    if (tree_nodes[tn].ends.size() > 1)
    {
      #define DIRECTED_DIFF // interpolates each point so it is exactly at max_distance, to avoid spurious splitting
        
      // 2. do floodfill on child roots to find if we have separate branches
      for (size_t i = 0; i<points.size(); i++)
        points[i].visited = false;   
      std::vector<int> new_ends;
      new_ends.push_back(tree_nodes[tn].ends[0]);
      for (size_t ijk = 0; ijk<new_ends.size(); ijk++)
      {
        int i = new_ends[ijk];
        points[i].visited = true;
        for (auto &j: tree_nodes[tn].ends)
        {
          if (points[j].visited)
            continue;
          #if defined DIRECTED_DIFF
          if (points[i].parent == -1 && points[j].parent == -1)
            std::cout << "something went wrong, end points should always have a parent" << std::endl;
          double dist1i = points[i].distance_to_ground;
          double dist0i = points[points[i].parent].distance_to_ground;
          double blendi = (max_dist_from_ground - dist0i) / (dist1i - dist0i);
          Eigen::Vector3d posi = points[points[i].parent].pos*(1.0-blendi) + points[i].pos*blendi;

          double dist1j = points[j].distance_to_ground;
          double dist0j = points[points[j].parent].distance_to_ground;
          double blendj = (max_dist_from_ground - dist0j) / (dist1j - dist0j);
          Eigen::Vector3d posj = points[points[j].parent].pos*(1.0-blendj) + points[j].pos*blendj;

          Eigen::Vector3d diff = posi - posj;
          #else
          Eigen::Vector3d diff = points[i].pos - points[j].pos;
          #endif
          if (diff.norm() < 2.0*tree_nodes[tn].radius)
          {
            new_ends.push_back(j);
            points[j].visited = true;
          }
        }        
      }

      if (new_ends.size() < tree_nodes[tn].ends.size()) // 3. if we are splitting then split and remove roots
      {
        extract_from_ends = true;
        nodes.clear(); // don't trust the found nodes as it is now two separate tree nodes

        TreesNode new_node = tree_nodes[tn];
        new_node.ends.clear();
        for (auto &node: tree_nodes[tn].ends)
        {
          if (std::find(new_ends.begin(), new_ends.end(), node) == new_ends.end())
            new_node.ends.push_back(node);
        }
        tree_nodes.push_back(new_node);
        if (verbose)
          std::cout << "connected ends: " << new_ends.size() << " so new node " << tree_nodes.size() << " with " << new_node.ends.size() << " ends" << std::endl;

        tree_nodes[tn].ends = new_ends;
      }
    }

    if (extract_from_ends) // 4. we have split the ends, so we need to extract the set of nodes in a backwards manner
    {
      for (auto &end: tree_nodes[tn].ends)
      {
        int node = points[end].parent;
        if (node == -1)
        {
          std::cout << "shouldn't be parentless here" << std::endl;
          continue;
        }
        while (node != -1)
        {
          if (std::find(nodes.begin(), nodes.end(), node) != nodes.end())
            break;
          nodes.push_back(node);
          if (std::find(tree_nodes[tn].roots.begin(), tree_nodes[tn].roots.end(), node) != nodes.end())
            break;
          node = points[node].parent;
        }
      }
      if (verbose)
        std::cout << "no roots, so working backwards has found " << nodes.size() << " nodes in total" << std::endl;
    }
    if (nodes.size() == 0)
    {
      if (verbose)
        std::cout << "bad node size zero" << std::endl;
      continue;
    }
    // Finally we have it, a set of nodes from which to get a centroid and radius estimation


    // find points in this segment, and the root points for the child segment
    tree_nodes[tn].centroid.setZero();
    for (auto &i: nodes)
      tree_nodes[tn].centroid += points[i].pos;

    if (nodes.size() > 0)
      tree_nodes[tn].centroid /= (double)nodes.size();
    
    Eigen::Vector3d dir(0,0,1);
    if (tree_nodes[tn].parent != -1)
      dir = (tree_nodes[tn].centroid - tree_nodes[tree_nodes[tn].parent].centroid).normalized();
    double rad = 0.0;
    for (auto &node: nodes)
    {
      Eigen::Vector3d offset = points[node].pos - tree_nodes[tn].centroid;
      Eigen::Vector3d p = offset - dir*offset.dot(dir);
      rad += p.squaredNorm();
    }
    if (nodes.size() > 4 || tree_nodes[tn].parent == -1)
    {
      rad /= (double)nodes.size();
      tree_nodes[tn].radius = std::sqrt(rad);
      std::cout << "estimated radius: " << tree_nodes[tn].radius << std::endl;
    }
    else
    {
      tree_nodes[tn].radius = tree_nodes[tree_nodes[tn].parent].radius;
    }

//    if (tree_nodes[tn].radius > 1.0)
//      tree_nodes[tn].radius = 0.05;//node_separation;

 //   tree_nodes[tn].radius = std::max(tree_nodes[tn].radius, 0.02);
    if (tree_nodes[tn].parent != -1)
      tree_nodes[tn].radius = std::min(tree_nodes[tn].radius, tree_nodes[tree_nodes[tn].parent].radius);

    // now add the single child for this particular tree node, assuming there are still ends
    if (tree_nodes[tn].ends.size() > 0)
    {
      TreesNode new_node;
      new_node.parent = (int)tn;
      new_node.roots = tree_nodes[tn].ends;
      new_node.radius = tree_nodes[tn].radius;
      new_node.min_dist_from_ground = max_dist_from_ground; // this is the only problem for 1-section-per-point... may affect something
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
        if ((tree_nodes[tree_node.parent].centroid - tree_node.centroid).norm() < 0.001)
          continue;
        if (tree_nodes[tree_node.parent].centroid.norm() < 0.1)
          continue;
        if (tree_node.centroid.norm() < 0.1)
          continue;
        starts.push_back(tree_nodes[tree_node.parent].centroid);
        ends.push_back(tree_node.centroid);
        radii.push_back(std::max(tree_node.radius, 0.01));
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
