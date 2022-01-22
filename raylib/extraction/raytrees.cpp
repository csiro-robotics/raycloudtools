// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include "../raydebugdraw.h"
#include "rayclusters.h"
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
  Vertex(const Eigen::Vector3d &pos) : pos(pos), edge_pos(0,0,0), parent(-1), distance_to_ground(inf), score(inf), visited(false) {}
  Eigen::Vector3d pos;
  Eigen::Vector3d edge_pos;
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
      if (i > 4000)
        break;
      if (points[i].parent != -1)
      {
        starts.push_back(points[points[i].parent].pos);
        ends.push_back(points[i].pos);
        Eigen::Vector3d col;
        col[0] = std::fmod(points[i].score, 1.0);
        col[1] = std::fmod(points[i].score/10.0,      1.0);
        col[2] = std::fmod(points[i].score/100.0, 1.0);
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
    if (!(sec%1000))
      std::cout << "sec " << sec << std::endl;
    // 1. Apply Leonardo's rule
    int par = sections[sec].parent;
    #if 0  // Leonardo's rule
    const double radius_change_scale = 1.1; // we're allowed to scale the total radius slightly each section
    if (par != -1 && sections[par].parent != -1) 
    {
      double rad = sections[sections[par].parent].radius;
      double child_rad = 0.0;
      for (auto &child: sections[sections[par].parent].children)
        child_rad += sqr(sections[child].radius);
      child_rad = std::sqrt(child_rad);
      rad = std::max(rad/radius_change_scale, std::min(child_rad, rad*radius_change_scale)); // try to head some way towards child_rad
      double scale = rad / child_rad;
      std::cout << "parent parent r " << rad << " has " << sections[sections[par].parent].children.size() << " children r " << child_rad << std::endl;
      for (auto &child: sections[sections[par].parent].children) // normalise children radii. Doesn't matter if we do this multiple times
        sections[child].radius *= scale;
    }
    #endif
    double thickness = 4.0*sections[sec].radius;
    if (par >= 0)
      thickness = 4.0*sections[par].radius;
    if (sec > 3000)
      break;

    double thickness_sqr = thickness*thickness;
    Eigen::Vector3d base(0,0,0);
    for (auto &root: sections[sec].roots)
      base += points[root].pos;
    base /= (double)sections[sec].roots.size();

    std::vector<int> nodes;
    std::cout << "tree " << sec << " radius used: " << thickness/4.0 << " roots: " << sections[sec].roots.size() << ", ends: " << sections[sec].ends.size() << std::endl;
    bool extract_from_ends = sections[sec].ends.size() > 0;
    if (!extract_from_ends)
    {
      nodes = sections[sec].roots;
      // 2. find all the points in this section:
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
    
      std::vector<int> all_ends = sections[sec].ends;
      // 3. do floodfill on child roots to find if we have separate branches
      for (auto &j: all_ends)
      {
        double dist1j = (points[j].pos - base).norm();
        double dist0j = (points[points[j].parent].pos - base).norm();
        double blendj = (thickness - dist0j) / (dist1j - dist0j);
        points[j].edge_pos = points[points[j].parent].pos*(1.0-blendj) + points[j].pos*blendj;
      }
      #define CLUSTER
      #if defined CLUSTER
      std::vector<Eigen::Vector3d> ps;
      std::vector<int> v_indices;
      for (auto &i: all_ends)
      {
        ps.push_back(points[i].edge_pos);
        v_indices.push_back(i);
      }
      std::vector< std::vector<int> > clusters = generateClusters(ps, 1.5*sections[sec].radius, 3.0*sections[sec].radius, true);
      for (auto &cluster: clusters)
      {
        for (auto &id: cluster)
          id = v_indices[id];
      }
      if (clusters.size() > 1)
      {
        std::cout << "first branch with " << clusters[0].size() << " / " << all_ends.size() << " points" << std::endl;
        extract_from_ends = true;
        nodes.clear(); // don't trust the found nodes as it is now two separate tree nodes
        sections[sec].ends = clusters[0]; // not quite, we need to translate back to real indices!
      }
      for (size_t i = 1; i<clusters.size(); i++)
      {
        std::cout << "subsequent branch with " << clusters[i].size() << " / " << all_ends.size() << " points" << std::endl;
        BranchSection new_node = sections[sec];
        new_node.ends = clusters[i];
        if (new_node.parent != -1)
          sections[new_node.parent].children.push_back((int)sections.size());
        sections.push_back(new_node);
      }
      #else
      int cc = -1;
      for (auto &end: all_ends)
      {
        cc++;
        #define DIRECTED_DIFF // interpolates each point so it is exactly at max_distance, to avoid spurious splitting
        if (points[end].visited)
          continue;        
        std::vector<int> new_ends;
        new_ends.push_back(end);
        for (size_t ijk = 0; ijk<new_ends.size(); ijk++)
        {
          int i = new_ends[ijk];
          points[i].visited = true;
          for (auto &j: all_ends)
          {
            if (points[j].visited)
              continue;
            Eigen::Vector3d diff = points[i].edge_pos - points[j].edge_pos;
            if (diff.norm() < 1.5*sections[sec].radius)
            {
              new_ends.push_back(j);
              points[j].visited = true;
            }
          }        
        }
        if (new_ends.size() < all_ends.size()) // 4. if we are splitting then split and remove roots
        {
          if (end == all_ends[0]) // if first clique
          {
            std::cout << "first branch with " << new_ends.size() << " / " << all_ends.size() << " points " << cc << std::endl;
            extract_from_ends = true;
            nodes.clear(); // don't trust the found nodes as it is now two separate tree nodes
            sections[sec].ends = new_ends;
          }
          else
          {
            std::cout << "subsequent branch with " << new_ends.size() << " / " << all_ends.size() << " points " << cc << std::endl;
            BranchSection new_node = sections[sec];
            new_node.ends = new_ends;
            if (new_node.parent != -1)
              sections[new_node.parent].children.push_back((int)sections.size());
            sections.push_back(new_node);
          }
        }
      }
      #endif
    }

    if (extract_from_ends) // 5. we have split the ends, so we need to extract the set of nodes in a backwards manner
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
    if (true)
    {
      std::vector<Eigen::Vector3d> debug_points;
      std::vector<double> shades;
      for (auto &i: sections[sec].ends)
      {
        debug_points.push_back(points[i].pos);
        shades.push_back(0.5);
      }      
      for (auto &i: nodes)
      {
        debug_points.push_back(points[i].pos);
        shades.push_back(1.0);
      }
      DebugDraw::instance()->drawCloud(debug_points, shades, 1);  
      std::vector<Eigen::Vector3d> starts;
      std::vector<Eigen::Vector3d> ends;
      for (auto &tree_node: sections)
      {
        if (sections[sec].tip == Eigen::Vector3d(0,0,0))
          continue;
        if (tree_node.parent >= 0)
        {
          if ((sections[tree_node.parent].tip - tree_node.tip).norm() < 0.001)
            continue;
          starts.push_back(sections[tree_node.parent].tip);
          ends.push_back(tree_node.tip);
        }
      }
      DebugDraw::instance()->drawLines(starts, ends);
    }

    // find points in this segment, and the root points for the child segment
    sections[sec].tip.setZero();
    for (auto &i: nodes)
      sections[sec].tip += points[i].pos;
    if (nodes.size() > 0)
      sections[sec].tip /= (double)nodes.size();

    
    // estimate radius
    if (par >= 0)
    {
      Eigen::Vector3d dir(0,0,1), prevdir(0,0,1);
      dir = (sections[sec].tip - sections[par].tip).normalized();
      if (sections[par].parent != -1)
        prevdir = (sections[par].tip - sections[sections[par].parent].tip).normalized();

      #define REAL_CENTROID
      #if defined REAL_CENTROID
      if (nodes.size() > 5)
      {
        Eigen::Vector3d mean_p(0,0,0);
        std::vector<Eigen::Vector3d> ps;
        Eigen::Vector3d vec(1,2,3);
        Eigen::Vector3d ax1 = dir.cross(vec).normalized();
        Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
        double n = 0;
        for (int ii = 0; ii<2; ii++)
        {
          std::vector<int> &node_list = ii==0 ? nodes : sections[sec].ends;
          for (auto &i: node_list) // one iteration of operation to find centre, using centroid, direction and radius estimation as a prior guess
          {
            Eigen::Vector3d pos = points[i].pos - sections[sec].tip;
            Eigen::Vector2d offset(ax1.dot(pos), ax2.dot(pos));
            Eigen::Vector2d xy = offset/sections[sec].radius;
            double l2 = xy.squaredNorm();
            Eigen::Vector3d point(xy[0], xy[1], 0.5*l2); // a paraboloid that has gradient 1 at 1
            ps.push_back(point);
            mean_p += point;   
            n++;  
          }
        }
        mean_p /= n; 
        struct Acc
        {
          Acc(){ x2 = y2 = xy = xz = yz = 0; }
          double x2, y2, xy, xz, yz;
        };
        Acc plane;
        for (auto &p: ps)
        {
          Eigen::Vector3d q = p - mean_p;
          plane.x2 += q[0]*q[0];
          plane.y2 += q[1]*q[1];
          plane.xy += q[0]*q[1];        
          plane.xz += q[0]*q[2];        
          plane.yz += q[1]*q[2];        
        }  
        const double eps = 1e-10;
        if (std::abs(plane.x2*plane.y2 - plane.xy*plane.xy) > eps && std::abs(plane.y2) > eps)
        {
          double A = (plane.xz*plane.y2 - plane.yz*plane.xy) / (plane.x2*plane.y2 - plane.xy*plane.xy);
          double B = (plane.yz - A * plane.xy) / plane.y2;

          Eigen::Vector2d shift(A,B);
          double shift2 = shift.squaredNorm();
          if (shift2 > 1.0) // don't shift more than one radius each iteration, for safety
            shift /= std::sqrt(shift2);

          std::cout << "shifting the centroid by " << shift.norm() * sections[sec].radius << " metres. New radius " << sections[sec].radius << std::endl;
          sections[sec].tip += (ax1*shift[0] + ax2*shift[1]) * sections[sec].radius;   
        }
      } 
      #endif   

      double rad = 0.0, rad_sqr = 0.0;
      for (auto &node: nodes)
      {
        Eigen::Vector3d offset = points[node].pos - sections[sec].tip;
        double dist_sqr = (offset - dir*offset.dot(dir)).squaredNorm();
        rad += std::sqrt(dist_sqr);
        rad_sqr += dist_sqr;
      }
      if (nodes.size() > 5)
      {
        double n = (double)nodes.size();
        double variance = (rad_sqr/n - sqr(rad / n)) * n/(n-4.0); // end part gives sample variance
        sections[sec].radius = std::sqrt(rad_sqr / n);
        std::cout << "estimated radius: " << sections[sec].radius << ", sigma: " << std::sqrt(variance) << std::endl;

        sections[sec].radius = std::max(sections[par].radius*0.75, sections[sec].radius);
        // constrain radius:
        if (std::sqrt(variance) > 0.025) 
        {
          sections[sec].radius = std::min(sections[sec].radius, sections[par].radius);
        }
      }
      else
      {
        sections[sec].radius = sections[par].radius;
        if (extract_from_ends) // multi-branching and not enough nodes for a reliable estimate... what do we do? use branch angle?
        {
          double sin_angle = dir.cross(prevdir).norm();
          std::cout << "not enough radius info on split, so using lateral change: " << sin_angle << " giving radius scale: " << 1.0/std::sqrt(1.0 + 2.0*sin_angle) << std::endl;
          sections[sec].radius /= std::sqrt(1.0 + 2.0*sin_angle); // increase the coefficient for more angle sensitivity 
        }
      }    
    }


    // now add the single child for this particular tree node, assuming there are still ends
    if (sections[sec].ends.size() > 0)
    {
      BranchSection new_node;
      new_node.parent = (int)sec;
      new_node.roots = sections[sec].ends;
      new_node.radius = sections[sec].radius;
      sections[sec].children.push_back((int)sections.size());
      sections.push_back(new_node);
    }    
  }
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<double> radii;
    for (auto &tree_node: sections)
    {
      if (tree_node.tip == Eigen::Vector3d(0,0,0))
        continue;
      if (tree_node.parent >= 0)
      {
        if ((sections[tree_node.parent].tip - tree_node.tip).norm() < 0.001)
          continue;
        starts.push_back(sections[tree_node.parent].tip);
        ends.push_back(tree_node.tip);
        radii.push_back(tree_node.radius);
      }
    }

    DebugDraw::instance()->drawLines(starts, ends);
    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
  }


  // generate children links
/*  for (unsigned int i = 0; i<sections.size(); i++)
  {
    if (sections[i].parent == -1)
      sections[sections[i].parent].children.push_back(i);
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
  ofs << "# trees file:" << std::endl;
  ofs << "x,y,z,radius,parent_id" << std::endl;
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
