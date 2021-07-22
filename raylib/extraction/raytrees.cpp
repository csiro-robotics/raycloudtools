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
  QueueNode(double distance_to_ground, double score, double radius, int index) : distance_to_ground(distance_to_ground), score(score), radius(radius), id(index) {}

  double distance_to_ground;
  double score;
  double radius;
  int id;
};

#define MINIMISE_SQUARE_DISTANCE // bad: end points are so distant that it creates separate branches
#define MINIMISE_ANGLE // works quite well in flowing along branches, but sometimes causes multi-branch problem, where radius was too small. 
class QueueNodeComparator 
{ 
public: 
    bool operator() (const QueueNode &p1, const QueueNode &p2) 
    { 
#if defined MINIMISE_SQUARE_DISTANCE || defined MINIMISE_ANGLE
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

Trees::Trees(const Cloud &cloud, const std::vector<std::pair<Eigen::Vector3d, double> > &trunks, bool verbose)
{
  std::vector<Vertex> points;  
  points.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
    if (cloud.rayBounded(i))
      points.push_back(Vertex(cloud.ends[i]));

	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;
  for (auto &trunk: trunks)
  {
    // find the lowest point in this trunk and put it a bit lower:
    double min_height = 1e10;
    for (auto &point: points)
    {
      Eigen::Vector3d dif = point.pos - trunk.first;
      dif[2] = 0.0;
      if (dif.squaredNorm() < 2.0*trunk.second*trunk.second)
        min_height = std::min(min_height, point.pos[2]);
    }
    Eigen::Vector3d root(trunk.first[0], trunk.first[1], min_height - trunk.second);
    BranchSection base;
    base.radius = trunk.second;
    closest_node.push(QueueNode(0, 0, base.radius, (int)points.size()));
    base.roots.push_back((int)points.size());
    points.push_back(Vertex(root));
    points.back().distance_to_ground = 0;
    points.back().score = 0;
    sections.push_back(base);
  }
  if (verbose)
  {
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }

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
  
  // 2b. climb up from lowest points...
	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!points[node.id].visited)
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        double dist = dists(i, node.id);
        double new_dist = node.distance_to_ground + dists(i, node.id)/node.radius;
        double new_score = 0;
        #if defined MINIMISE_SQUARE_DISTANCE
        dist *= dist;
        #endif
        #if defined MINIMISE_ANGLE
        Eigen::Vector3d dif = (points[child].pos - points[node.id].pos).normalized();
        Eigen::Vector3d dir(0,0,1);
        if (points[node.id].parent != -1)
          dir = (points[node.id].pos - points[points[node.id].parent].pos).normalized();
        const double power = 4.0;
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
					closest_node.push(QueueNode(points[child].distance_to_ground, points[child].score, node.radius, child));
				}
			}
		  points[node.id].visited = true;
		}
	}

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
        int slot = (int)(points[i].distance_to_ground / 0.2);
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


  // now trace from root tree nodes upwards, getting node centroids
  // a tree_node is a segment, and these are added as we iterate through the list
  for (size_t sec = 0; sec < sections.size(); sec++)
  {   
    double thickness = 4.0*sections[sec].radius;
    if (sections[sec].parent >= 0)
      thickness = 4.0*sections[sections[sec].parent].radius;
    double thickness_sqr = thickness*thickness;
    Eigen::Vector3d base(0,0,0);
    for (auto &root: sections[sec].roots)
      base += points[root].pos;
    base /= (double)sections[sec].roots.size();

    std::vector<int> nodes;
    std::cout << "tree " << sec << " roots: " << sections[sec].roots.size() << ", ends: " << sections[sec].ends.size() << std::endl;
    bool extract_from_ends = sections[sec].ends.size() > 0;
    if (!extract_from_ends)
    {
      nodes = sections[sec].roots;
      // 1. find all the points in this tree node:
      for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
      {
        int i = nodes[ijk];
        for (auto &child: children[i])
        {
          if ((points[child].pos - base).squaredNorm() < thickness_sqr) // in same slot, so accumulate
            nodes.push_back(child); // so we recurse on this child too
          else 
          {
            sections[sec].ends.push_back(child); 
          }
        }
      }
      if (verbose)
        std::cout << "no ends, so found " << sections[sec].ends.size() << " ends" << std::endl;
    }
    // TODO: what if we have found 0 ends here? i.e. the end of the branch...?
    if (sections[sec].ends.size() > 1)
    {
      #define DIRECTED_DIFF // interpolates each point so it is exactly at max_distance, to avoid spurious splitting
        
      // 2. do floodfill on child roots to find if we have separate branches
      for (size_t i = 0; i<points.size(); i++)
        points[i].visited = false;   
      std::vector<int> new_ends;
      new_ends.push_back(sections[sec].ends[0]);
      for (size_t ijk = 0; ijk<new_ends.size(); ijk++)
      {
        int i = new_ends[ijk];
        points[i].visited = true;
        for (auto &j: sections[sec].ends)
        {
          if (points[j].visited)
            continue;
          #if defined DIRECTED_DIFF
          if (points[i].parent == -1 && points[j].parent == -1)
            std::cout << "something went wrong, end points should always have a parent" << std::endl;
          double dist1i = (points[i].pos - base).norm();
          double dist0i = (points[points[i].parent].pos - base).norm();
          double blendi = (thickness - dist0i) / (dist1i - dist0i);
          Eigen::Vector3d posi = points[points[i].parent].pos*(1.0-blendi) + points[i].pos*blendi;

          double dist1j = (points[j].pos - base).norm();
          double dist0j = (points[points[j].parent].pos - base).norm();
          double blendj = (thickness - dist0j) / (dist1j - dist0j);
          Eigen::Vector3d posj = points[points[j].parent].pos*(1.0-blendj) + points[j].pos*blendj;

          Eigen::Vector3d diff = posi - posj;
          #else
          Eigen::Vector3d diff = points[i].pos - points[j].pos;
          #endif
          if (diff.norm() < 3.0*sections[sec].radius)
          {
            new_ends.push_back(j);
            points[j].visited = true;
          }
        }        
      }

      if (new_ends.size() < sections[sec].ends.size()) // 3. if we are splitting then split and remove roots
      {
        extract_from_ends = true;
        nodes.clear(); // don't trust the found nodes as it is now two separate tree nodes

        BranchSection new_node = sections[sec];
        new_node.ends.clear();
        for (auto &node: sections[sec].ends)
        {
          if (std::find(new_ends.begin(), new_ends.end(), node) == new_ends.end())
            new_node.ends.push_back(node);
        }
        if (verbose)
          std::cout << "connected ends: " << new_ends.size() << " so new node " << sections.size() << " with " << new_node.ends.size() << " ends" << std::endl;
        
        sections.push_back(new_node);
        sections[sec].ends = new_ends;
      }
    }
    if (extract_from_ends) // 4. we have split the ends, so we need to extract the set of nodes in a backwards manner
    {
      for (auto &end: sections[sec].ends)
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
          if (std::find(sections[sec].roots.begin(), sections[sec].roots.end(), node) != sections[sec].roots.end())
            break;
          node = points[node].parent;
        }
      }
      if (verbose)
        std::cout << "extract from ends, so working backwards has found " << nodes.size() << " nodes in total" << std::endl;
    }
    // Finally we have it, a set of nodes from which to get a centroid and radius estimation

    // find points in this segment, and the root points for the child segment
    Eigen::Vector3d centroid(0,0,0);
    for (auto &i: nodes)
      centroid += points[i].pos;
    if (nodes.size() > 0)
      centroid /= (double)nodes.size();
    sections[sec].tip = centroid; 

    Eigen::Vector3d dir(0,0,1);
    if (sections[sec].parent != -1)
      dir = (sections[sec].tip - sections[sections[sec].parent].tip).normalized();
    
    // use the section centroid for estimating the radius
    double rad = 0.0;
    for (auto &node: nodes)
    {
      Eigen::Vector3d offset = points[node].pos - centroid;
      Eigen::Vector3d p = offset - dir*offset.dot(dir);
      rad += p.squaredNorm();
    }
    if (nodes.size() > 5 || sections[sec].parent == -1)
    {
      rad /= (double)nodes.size();
      sections[sec].radius = std::sqrt(rad);
      std::cout << "estimated radius: " << sections[sec].radius << std::endl;
    }
    else
    {
      sections[sec].radius = sections[sections[sec].parent].radius;
    }

    if (sections[sec].parent != -1)
    {
      sections[sec].radius = std::min(sections[sec].radius, sections[sections[sec].parent].radius);
      sections[sec].radius = std::max(sections[sec].radius, 0.7*sections[sections[sec].parent].radius);
    }

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections[sec].ends.size() > 0)
    {
      BranchSection new_node;
      new_node.parent = (int)sec;
      new_node.roots = sections[sec].ends;
      new_node.radius = sections[sec].radius;
      sections.push_back(new_node);
    }    
  }
  // remove node 0
  for (auto &tree_node: sections)
  {
    if (tree_node.parent >= 0 && sections[tree_node.parent].parent == -1)
      tree_node.parent = -2;
  }
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<double> radii;
    for (auto &tree_node: sections)
    {
      if (tree_node.parent >= 0)
      {
        if ((sections[tree_node.parent].tip - tree_node.tip).norm() < 0.001)
          continue;
        if (sections[tree_node.parent].tip.norm() < 0.1)
          continue;
        if (tree_node.tip.norm() < 0.1)
          continue;
        starts.push_back(sections[tree_node.parent].tip);
        ends.push_back(tree_node.tip);
        radii.push_back(std::max(tree_node.radius, 0.01));
      }
    }
    DebugDraw::instance()->drawLines(starts, ends);
    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
  }


  // generate children links
  for (unsigned int i = 0; i<sections.size(); i++)
  {
    if (sections[i].parent >= 0)
    {
      sections[sections[i].parent].children.push_back(i);
    }
    else if (sections[i].parent == -2)
      root_nodes.push_back(i);
  }
  // generate local parent links
  for (auto &root: root_nodes)
  {
    int child_id = 0;
    sections[root].id = child_id++;
    std::vector<int> children = sections[root].children;
    for (unsigned int c = 0; c<children.size(); c++)
    {
      sections[children[c]].id = child_id++;
      for (auto i: sections[children[c]].children)
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
    ofs << sections[root].tip[0] << "," << sections[root].tip[1] << "," << sections[root].tip[2] << "," << sections[root].radius << ",-1";

    std::vector<int> children = sections[root].children;
    for (unsigned int c = 0; c<children.size(); c++)
    {
      BranchSection &node = sections[children[c]];
      ofs << ", " << node.tip[0] << "," << node.tip[1] << "," << node.tip[2] << "," << node.radius << "," << sections[node.parent].id;
      for (auto i: sections[children[c]].children)
        children.push_back(i);
    }
    ofs << std::endl;
  }
  return true;
}

} // namespace ray
