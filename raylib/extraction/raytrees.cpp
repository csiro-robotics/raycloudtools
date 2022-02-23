// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include "../raydebugdraw.h"
#include "rayclusters.h"
#include <nabo/nabo.h>

namespace ray
{
double radFromLength(double length, const TreesParams &params)
{
  if (length / params.linear_range)
    return length / params.length_to_radius;
  double radius_at_min = params.linear_range / params.length_to_radius;
  double unscaled_rad = std::pow(radius_at_min, params.radius_exponent);
  // length = scale_factor * rad^radius_exponent
  return std::pow(length * unscaled_rad / params.linear_range, 1.0/params.radius_exponent);
}

Trees::Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose)
{
  std::vector<Vertex> points;  
  std::vector< std::vector<int> > roots_list = getRootsAndSegment(points, cloud, mesh, params.max_diameter, params.distance_limit, params.height_min, params.gravity_factor);
  const int orig_points = (int)points.size() - (int)mesh.vertices().size();
  if (verbose)
    cloud.save("test_output.ply");

  for (size_t i = 0; i<roots_list.size(); i++)
  {
    BranchSection base;
    base.roots = roots_list[i];
    for (auto &root: base.roots)
      if (root < orig_points)
        std::cout << "error: " << root << std::endl;
    sections.push_back(base);
  }
  if (verbose)
  {
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }

  // now get distance to end
  for (size_t i = 0; i<points.size(); i++)
  {
    int id = (int)i;
    int parent = points[i].parent;
    while (parent != -1)
    {
      double end_dist = points[id].distance_to_end + (points[id].distance_to_ground - points[parent].distance_to_ground);
      if (points[parent].distance_to_end > end_dist)
        break;
      points[parent].distance_to_end = end_dist;
      id = parent;
      parent = points[id].parent;
    }
  }

  for (int i = (int)sections.size()-1; i>=0; i--)
  {
    for (auto &root: sections[i].roots)
      sections[i].max_distance_to_end = std::max(sections[i].max_distance_to_end, points[root].distance_to_end);
    sections[i].radius = radFromLength(sections[i].max_distance_to_end, params);
 //   std::cout << "initial end distances: " << sections[i].roots.size() << ": " << sections[i].max_distance_to_end << " rad: " << sections[i].radius << std::endl;
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
  std::vector<int> section_ids(points.size(), -1);

  // now trace from root tree nodes upwards, getting node centroids
  // a tree_node is a segment, and these are added as we iterate through the list
  for (size_t sec = 0; sec < sections.size(); sec++)
  {   
    if (!(sec%10000))
      std::cout << "generating segment " << sec << std::endl;
    int par = sections[sec].parent;

    double max_radius = radFromLength(sections[sec].max_distance_to_end, params);
    if (max_radius < params.minimum_radius)
    {
      sections[sec].tip.setZero();
      for (auto &i: sections[sec].ends)
        sections[sec].tip += points[i].pos;
      size_t sum = sections[sec].ends.size();
      if (sections[sec].ends.empty())
      {
        for (auto &i: sections[sec].roots)
          sections[sec].tip += points[i].pos;
        sum = sections[sec].roots.size();
      }
      if (sum > 0)
        sections[sec].tip /= (double)sum;
      continue;
    }

    double thickness = params.cylinder_length_to_width * max_radius;
    double thickness_sqr = thickness*thickness;
    std::vector<int> nodes;
    bool extract_from_ends = sections[sec].ends.size() > 0;
    if (!extract_from_ends)
    {
      Eigen::Vector3d base(0,0,0);
      for (auto &root: sections[sec].roots)
        base += points[root].pos;
      base /= (double)sections[sec].roots.size();

      nodes = sections[sec].roots;
      // 2. find all the points (and the end points) for this section:
      for (unsigned int ijk = 0; ijk<nodes.size(); ijk++)
      {
        int i = nodes[ijk];
        for (auto &child: children[i])
        {
          double dist_sqr = par == -1 ? sqr(points[child].pos[2]-base[2]) : (points[child].pos - base).squaredNorm();
          if (dist_sqr < thickness_sqr) // in same slot, so accumulate
            nodes.push_back(child); // so we recurse on this child too
          else 
            sections[sec].ends.push_back(child); 
        }
      }
    
      std::vector<int> all_ends = sections[sec].ends;
      // 3. cluster child roots to find if we have separate branches
      // first, get interpolated edge points
      for (auto &j: all_ends)
      {
        double dist1j = par==-1 ? points[j].pos[2]-base[2] : (points[j].pos-base).norm();
        double dist0j = par==-1 ? points[points[j].parent].pos[2]-base[2] : (points[points[j].parent].pos-base).norm();
        double blendj = (thickness - dist0j) / (dist1j - dist0j);
        points[j].edge_pos = points[points[j].parent].pos*(1.0-blendj) + points[j].pos*blendj;
      }
      std::vector<Eigen::Vector3d> ps;
      std::vector<int> v_indices;
      for (auto &i: all_ends)
      {
        ps.push_back(points[i].edge_pos);
        v_indices.push_back(i);
      }
      std::vector< std::vector<int> > clusters = generateClusters(ps, params.gap_ratio*max_radius, params.span_ratio*max_radius, true);
      for (auto &cluster: clusters)
      {
        for (auto &id: cluster)
          id = v_indices[id];
      }
      size_t min_size = 1;
      if (par == -1)
      {
        // remove children that are too small
        for (int i = (int)clusters.size()-1; i>=0; i--)
        {
          double max_dist = 0.0;
          for (auto &end: clusters[i])
          {
            max_dist = std::max(max_dist, points[end].distance_to_end);
          }
          if (max_dist < params.height_min)
          {
            clusters[i] = clusters.back();
            clusters.pop_back();
            min_size = 0;            
          }
        }
      }
      std::vector<double> max_distances(clusters.size());
      double maxmax = -1;
      int maxi = -1;
      for (size_t i = 0; i<clusters.size(); i++)
      {
        max_distances[i] = 0;
        for (auto &end: clusters[i])
          max_distances[i] = std::max(max_distances[i], points[end].distance_to_end);
        if (max_distances[i] > maxmax)
        {
          maxmax = max_distances[i];
          maxi = (int)i;
        }
      }
      if (clusters.size() > min_size)
      {
        if (maxi == -1)
          std::cout << "error: bad maxi" << std::endl;
        extract_from_ends = true;
        nodes.clear(); // don't trust the found nodes as it is now two separate tree nodes
        sections[sec].ends = clusters[maxi]; 
        sections[sec].max_distance_to_end = max_distances[maxi] + thickness;
        max_radius = radFromLength(sections[sec].max_distance_to_end, params);
      }
      for (size_t i = 0; i<clusters.size(); i++)
      {
        if ((int)i == maxi)
          continue;
        BranchSection new_node = sections[sec];
        new_node.max_distance_to_end = max_distances[i] + thickness;
        double maxrad = radFromLength(new_node.max_distance_to_end, params);
        if (maxrad > params.minimum_radius)
        {
          new_node.ends = clusters[i];      
          if (par != -1)
            sections[par].children.push_back((int)sections.size());
          sections.push_back(new_node);
        }
      }
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
        if (par == -1)
        {
          section_ids[nodes.back()] = (int)sec;
        }
      }
    }
    else if (par == -1)
    {
      for (auto &root: sections[sec].roots)
      {
        section_ids[root] = (int)sec;
      }    
    }

    // find points in this segment, and the root points for the child segment
    sections[sec].tip.setZero();
    if (nodes.empty())
      std::cout << "this shouldn't happen!!!! empty nodes" << std::endl;
    auto &list = nodes.empty() ? (sections[sec].ends.empty() ? sections[sec].roots : sections[sec].ends) : nodes;
    for (auto &i: list)
      sections[sec].tip += points[i].pos;
    if (list.size() > 0)
      sections[sec].tip /= (double)list.size();

    // estimate radius
    {
      Eigen::Vector3d dir(0,0,1);
      if (par >= 0)
        dir = (sections[sec].tip - sections[par].tip).normalized();

      #define REAL_CENTROID
      #if defined REAL_CENTROID
      Eigen::Vector3d mean_p(0,0,0);
      std::vector<Eigen::Vector3d> ps;
      Eigen::Vector3d vec(1,2,3);
      Eigen::Vector3d ax1 = dir.cross(vec).normalized();
      Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
      double n = 0;
      int start_ii = par >= 0 || sections[sec].ends.empty() ? 0 : 1;
      for (int ii = start_ii; ii<2; ii++)
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
      if (n > 5)
      {
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
          if (par >= 0 && shift2 > 1.0) // don't shift more than one radius each iteration, for safety
            shift /= std::sqrt(shift2);

          sections[sec].tip += (ax1*shift[0] + ax2*shift[1]) * sections[sec].radius;   
        }
      } 
      else if (n > 0)
      {
        sections[sec].tip += (ax1*mean_p[0] + ax2*mean_p[1]) * sections[sec].radius; 
      }
      #endif   

      if (par == -1)
      {
        double n = 0, rad = 0;
        const auto &list = sections[sec].ends.empty() ? nodes : sections[sec].ends;
        for (auto &id: list)
        {
          Eigen::Vector3d offset = points[id].pos - sections[sec].tip;
          rad += (offset - dir*offset.dot(dir)).norm();
          n++;
        }
        double radius = list.size() < 2 ? max_radius : rad / n;
        sections[sec].radius = std::min(radius, max_radius);
      }
      else
      {
        double n = 4.0;
        double rad = n * sections[par].radius;
        for (auto &node: nodes)
        {
          Eigen::Vector3d offset = points[node].pos - sections[sec].tip;
          rad += (offset - dir*offset.dot(dir)).norm();
          n++;
        }
        sections[sec].radius = std::min( std::min(rad / n, max_radius), sections[par].radius);
      }
    }

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections[sec].ends.size() > 0)
    {
      BranchSection new_node;
      new_node.parent = (int)sec;
      new_node.roots = sections[sec].ends;
      new_node.max_distance_to_end = 0.0;
      for (auto &root: new_node.roots)
        new_node.max_distance_to_end = std::max(new_node.max_distance_to_end, points[root].distance_to_end);
      max_radius = radFromLength(new_node.max_distance_to_end, params);
      new_node.radius = std::min(sections[sec].radius, max_radius);
      if (new_node.radius > params.minimum_radius) // if it is the first node, then we need a second noded
      {
        sections[sec].children.push_back((int)sections.size());

        new_node.tip.setZero();
        for (auto &end: new_node.roots)
          new_node.tip += points[end].pos;
        if (new_node.roots.size() > 0)
          new_node.tip /= (double)new_node.roots.size();

        sections.push_back(new_node);
      }
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
    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
  }
  // generate local parent links
  int num = 0;
  for (auto &section: sections)
  {
    if (section.parent >= 0)
      continue;
    num++;
    int child_id = 0;
    section.id = child_id++;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c<children.size(); c++)
    {
      sections[children[c]].id = child_id++;
      for (auto i: sections[children[c]].children)
        children.push_back(i);
    }
  }
  std::cout << num << " trees saved" << std::endl;

  int j = -1;
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    RGBA &colour = cloud.colours[i];
    if (cloud.rayBounded(i))
    {
      j++;
      int root = points[j].root;
      if (root == -1)
      {
        colour.red = colour.green = colour.blue = 0;
        continue;
      }
      int seg = section_ids[root];
      if (seg == -1)
      {
        colour.red = colour.green = colour.blue = 0;
        continue;
      }
      srand(1 + seg);
      colour.red   = uint8_t(50 + rand()%205);
      colour.green = uint8_t(50 + rand()%205);
      colour.blue  = uint8_t(50 + rand()%205);
    }
    else
    {
      colour.red = colour.green = colour.blue = 0;
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
  ofs << "# tree structure file:" << std::endl;
  ofs << "x,y,z,radius,parent_id" << std::endl;
  for (const auto &section: sections)
  {
    if (section.parent >= 0 || section.children.empty())
      continue;
    ofs << section.tip[0] << "," << section.tip[1] << "," << section.tip[2] << "," << section.radius << ",-1";

    std::vector<int> children = section.children;
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
