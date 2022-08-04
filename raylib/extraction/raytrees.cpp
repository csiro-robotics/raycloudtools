// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include <nabo/nabo.h>
#include "../raydebugdraw.h"
#include "rayclusters.h"

namespace ray
{
/// calculate a branch radius from its length. We use allometric scaling based on a radius exponent and a 
/// linear range for the final taper of the branch
/// See "Allometric patterns in Acer platanoides (Aceraceae) branches" for real world data
double radFromLength(double length, const TreesParams &params)
{
  if (length / params.linear_range)
  {
    return length / params.length_to_radius;
  }
  double radius_at_min = params.linear_range / params.length_to_radius;
  double unscaled_rad = std::pow(radius_at_min, params.radius_exponent);
  // length = scale_factor * rad^radius_exponent
  return std::pow(length * unscaled_rad / params.linear_range, 1.0 / params.radius_exponent);
}

/// The main reconstruction algorithm
/// It is based on finding the shortest paths using Djikstra's algorithm, followed
/// by an agglomeration of paths, with repeated splitting from root to tips
Trees::Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose)
{
  // firstly, get the full set of shortest paths from ground to tips, and the set of roots
  std::vector<Vertex> points;
  std::vector<std::vector<int>> roots_list = getRootsAndSegment(
    points, cloud, mesh, params.max_diameter, params.distance_limit, params.height_min, params.gravity_factor);
  const int orig_points = (int)points.size() - (int)mesh.vertices().size();

  // Now we want to convert these paths into a set of branch sections, from root to tips
  // splitting as we go up. 

  // first we initialise to one branch section per root
  std::vector<int> section_ids(points.size(), -1);
  for (size_t i = 0; i < roots_list.size(); i++)
  {
    BranchSection base;
    base.roots = roots_list[i];
    for (auto &root : base.roots)
    {
      section_ids[root] = (int)sections.size();
      if (root < orig_points)
        std::cout << "error: " << root << std::endl;
    }
    sections.push_back(base);
  }
  if (verbose)
  {
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }

  // now calculate the maximum distance to tip for every point in the cloud
  // this will be helpful in the next step
  for (size_t i = 0; i < points.size(); i++)
  {
    int id = (int)i;
    int parent = points[i].parent;
    while (parent != -1)
    {
      double end_dist =
        points[id].distance_to_end + (points[id].distance_to_ground - points[parent].distance_to_ground);
      if (points[parent].distance_to_end > end_dist)
        break;
      points[parent].distance_to_end = end_dist;
      id = parent;
      parent = points[id].parent;
    }
  }

  // we can now find the maximum distance to tip per root section, which allows us to estimate the 
  // section radius for all the root sections. 
  for (int i = (int)sections.size() - 1; i >= 0; i--)
  {
    for (auto &root : sections[i].roots)
      sections[i].max_distance_to_end = std::max(sections[i].max_distance_to_end, points[root].distance_to_end);
    sections[i].radius = radFromLength(sections[i].max_distance_to_end, params);
  }

  // the set of child points will be useful later, so generate them now
  std::vector<std::vector<int>> children(points.size());
  for (size_t i = 0; i < points.size(); i++)
  {
    if (points[i].parent != -1)
      children[points[i].parent].push_back((int)i);
  }
  for (size_t i = 0; i < points.size(); i++)
  {
    points[i].visited = false;
  }

  // now trace from root tree nodes upwards, getting node centroids
  // create new BranchSections as we go
  for (size_t sec = 0; sec < sections.size(); sec++)
  {
    if (!(sec % 10000))
      std::cout << "generating segment " << sec << std::endl;
    int par = sections[sec].parent;

    // any estimate of branch radius is bounded from above by an estimate from the branch's length
    // estimates from length are more reliable on small branches, and less reliable on large, well
    // observed trunks. 
    double max_radius = radFromLength(sections[sec].max_distance_to_end, params);

    // this branch can come to an end as it is now too small
    if (max_radius < 0.5 * params.min_diameter)
    {
      // before we finish, we need to calculate the tip of this branch, from its end points
      sections[sec].tip.setZero();
      for (auto &i : sections[sec].ends) 
      {
        sections[sec].tip += points[i].pos;
      }
      size_t sum = sections[sec].ends.size();
      // if it doesnt have end points then we'll have to use the root points
      if (sections[sec].ends.empty()) 
      {
        for (auto &i : sections[sec].roots) 
        {
          sections[sec].tip += points[i].pos;
        }
        sum = sections[sec].roots.size();
      }
      if (sum > 0)
      {
        sections[sec].tip /= (double)sum;
      }
      continue;
    }

    double thickness = params.cylinder_length_to_width * max_radius;
    double thickness_sqr = thickness * thickness;
    std::vector<int> nodes; // all the points in the section
    bool extract_from_ends = sections[sec].ends.size() > 0;
    // if the branch section has no end points recorded, then we need to examine this branch to 
    // find end points and potentially branch (bifurcate)  
    if (!extract_from_ends)
    {
      // get a base position for the section
      Eigen::Vector3d base(0, 0, 0);
      for (auto &root : sections[sec].roots) 
      {
        base += points[root].pos;
      }
      base /= (double)sections[sec].roots.size();

      nodes = sections[sec].roots;
      // 2. find all the points (and the end points) for this section:
      for (unsigned int ijk = 0; ijk < nodes.size(); ijk++)
      {
        int i = nodes[ijk];
        for (auto &child : children[i])
        {
          double dist_sqr = par == -1 ? sqr(points[child].pos[2] - base[2]) : (points[child].pos - base).squaredNorm();
          if (dist_sqr < thickness_sqr)  // in same slot, so accumulate
          {
            nodes.push_back(child);      // so we recurse on this child too
          }
          else
          {
            sections[sec].ends.push_back(child);
          }
        }
      }

      std::vector<int> all_ends = sections[sec].ends;
      // 3. cluster end points to find if we have separate branches
      // first, get interpolated edge points. i.e. interpolate between two connected points inside and outside
      // the branch section's cylinder, so that the set edge_pos are right on the top boundary of the section
      for (auto &j : all_ends)
      {
        double dist1j = par == -1 ? points[j].pos[2] - base[2] : (points[j].pos - base).norm();
        double dist0j =
          par == -1 ? points[points[j].parent].pos[2] - base[2] : (points[points[j].parent].pos - base).norm();
        double blendj = (thickness - dist0j) / (dist1j - dist0j);
        points[j].edge_pos = points[points[j].parent].pos * (1.0 - blendj) + points[j].pos * blendj;
      }
      // convert to a structure that is better for the cluster function
      std::vector<Eigen::Vector3d> ps;
      std::vector<int> v_indices;
      for (auto &i : all_ends)
      {
        ps.push_back(points[i].edge_pos);
        v_indices.push_back(i);
      }
      // cluster these end points based on two separation criteria (gap_ratio and span_ratio)
      std::vector<std::vector<int>> clusters =
        generateClusters(ps, params.gap_ratio * max_radius, params.span_ratio * max_radius);
      // adjust back to global ids
      for (auto &cluster : clusters)
      {
        for (auto &id : cluster) 
        {
          id = v_indices[id];
        }
      }

      size_t min_size = 1;
      if (par == -1) // if this is the root section
      {
        // then remove children that are smaller than the minimum tree height
        for (int i = (int)clusters.size() - 1; i >= 0; i--)
        {
          double max_dist = 0.0;
          for (auto &end : clusters[i])
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

      // find the maximum distance (to tip) for each cluster
      std::vector<double> max_distances(clusters.size());
      double maxmax = -1;
      int maxi = -1;
      for (size_t i = 0; i < clusters.size(); i++)
      {
        max_distances[i] = 0;
        for (auto &end : clusters[i]) 
        {
          max_distances[i] = std::max(max_distances[i], points[end].distance_to_end);
        }
        if (max_distances[i] > maxmax)
        {
          maxmax = max_distances[i];
          maxi = (int)i;
        }
      }

      // if we are bifurcating then set the current branch section to the cluster with the
      // maximum distance to tip (i.e. the longest branch). 
      if (clusters.size() > min_size) 
      {
        if (maxi == -1)
          std::cout << "error: bad maxi" << std::endl;
        extract_from_ends = true;
        nodes.clear();  // don't trust the found nodes as it is now two separate tree nodes
        sections[sec].ends = clusters[maxi];
        sections[sec].max_distance_to_end = max_distances[maxi] + thickness;
        max_radius = radFromLength(sections[sec].max_distance_to_end, params);
      }
      // for all other clusters, add new sections to the list...
      for (size_t i = 0; i < clusters.size(); i++)
      {
        if ((int)i == maxi)
          continue;
        BranchSection new_node = sections[sec];
        new_node.max_distance_to_end = max_distances[i] + thickness;
        double maxrad = radFromLength(new_node.max_distance_to_end, params);
        if (maxrad > 0.5 * params.min_diameter) // but only add if they are large enough
        {
          // we only specify the end points at this stage. They will therefore enter the 
          // extract_from_ends block below when the sec for loop gets to their section
          new_node.ends = clusters[i]; 
          if (par != -1)
          {
            sections[par].children.push_back((int)sections.size());
          }
          sections.push_back(new_node);
        }
      }
    }

    // we have split the ends, so we need to extract the set of nodes in a backwards manner
    if (extract_from_ends)  
    {
      for (auto &end : sections[sec].ends)
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
          {
            break;
          }
          nodes.push_back(node); // fill in the nodes in this branch section
          if (std::find(sections[sec].roots.begin(), sections[sec].roots.end(), node) != sections[sec].roots.end())
          {
            break;
          }
          node = points[node].parent;
        }
      }
    }

    if (nodes.empty())
    {
      std::cout << "error: there shouldn't be empty nodes at this point in the processing" << std::endl;
    }
    // get the points in this segment
    auto &list = nodes.empty() ? (sections[sec].ends.empty() ? sections[sec].roots : sections[sec].ends) : nodes;

    // use the nodes to estimate a tip location (which is the mean of the node points)
    sections[sec].tip.setZero();
    for (auto &i : list) 
    {
      sections[sec].tip += points[i].pos;
    }
    if (list.size() > 0)
    {
      sections[sec].tip /= (double)list.size();
    }

    // estimate branch section radius
    {
      // get section's direction
      Eigen::Vector3d dir(0, 0, 1);
      if (par >= 0)
      {
        dir = (sections[sec].tip - sections[par].tip).normalized();
      }

#define REAL_CENTROID // finds a new centroid that is robust to branches scanned from a single side
#if defined REAL_CENTROID
      Eigen::Vector3d mean_p(0, 0, 0);
      std::vector<Eigen::Vector3d> ps;
      Eigen::Vector3d vec(1, 2, 3);
      // obtain two orthogonal planes to the section's direction vector
      Eigen::Vector3d ax1 = dir.cross(vec).normalized();
      Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
      double n = 0;
      // iterate over the nodes and the ends points...
      int start_ii = par >= 0 || sections[sec].ends.empty() ? 0 : 1;
      for (int ii = start_ii; ii < 2; ii++)
      {
        std::vector<int> &node_list = ii == 0 ? nodes : sections[sec].ends;
        // project points into a local space parabola
        for (auto &i : node_list)  
        {
          Eigen::Vector3d pos = points[i].pos - sections[sec].tip;
          Eigen::Vector2d offset(ax1.dot(pos), ax2.dot(pos));
          Eigen::Vector2d xy = offset / sections[sec].radius;
          double l2 = xy.squaredNorm();
          Eigen::Vector3d point(xy[0], xy[1], 0.5 * l2);  // a paraboloid that has gradient 1 at 1
          ps.push_back(point);
          mean_p += point;
          n++;
        }
      }
      mean_p /= n;
      if (n > 5) // assuming there are sufficient points for a resonable guess
      {
        // accumulation structure for plane least squares fitting
        struct Acc
        {
          Acc() { x2 = y2 = xy = xz = yz = 0; }
          double x2, y2, xy, xz, yz;
        };
        Acc plane;
        for (auto &p : ps)
        {
          // fill in the parameters to estimate the plane of best fit to the paraboloid
          Eigen::Vector3d q = p - mean_p;
          plane.x2 += q[0] * q[0];
          plane.y2 += q[1] * q[1];
          plane.xy += q[0] * q[1];
          plane.xz += q[0] * q[2];
          plane.yz += q[1] * q[2];
        }
        const double eps = 1e-10;
        // is the plane is well-determined
        if (std::abs(plane.x2 * plane.y2 - plane.xy * plane.xy) > eps && std::abs(plane.y2) > eps)
        {
          // extract the local plane parameters
          double A = (plane.xz * plane.y2 - plane.yz * plane.xy) / (plane.x2 * plane.y2 - plane.xy * plane.xy);
          double B = (plane.yz - A * plane.xy) / plane.y2;

          Eigen::Vector2d shift(A, B);
          double shift2 = shift.squaredNorm();
          if (par >= 0 && shift2 > 1.0)  // don't shift more than one radius each iteration, for safety
          { 
            shift /= std::sqrt(shift2);
          }

          // apply the plane parameters as a world-space shift in the branch section tip position
          sections[sec].tip += (ax1 * shift[0] + ax2 * shift[1]) * sections[sec].radius;
        }
      }
      else if (n > 0) // if only a few points then use the mean instead
      {
        sections[sec].tip += (ax1 * mean_p[0] + ax2 * mean_p[1]) * sections[sec].radius;
      }
#endif

      // now find the segment radius
      if (par == -1) // if this is the root segment
      {
        double n = 0, rad = 0;
        // then get the mean radius
        const auto &list = sections[sec].ends.empty() ? nodes : sections[sec].ends;
        for (auto &id : list)
        {
          Eigen::Vector3d offset = points[id].pos - sections[sec].tip;
          rad += (offset - dir * offset.dot(dir)).norm();
          n++;
        }
        double radius = list.size() < 2 ? max_radius : rad / n;
        sections[sec].radius = std::min(radius, max_radius);
      }
      else // for non-root segments 
      {
        // use the parent radius as a prior with a weight of 4 points
        // this avoids spurious radius estimations when the number of points is
        // small
        double n = 4.0; 
        double rad = n * sections[par].radius;
        for (auto &node : nodes)
        {
          Eigen::Vector3d offset = points[node].pos - sections[sec].tip;
          rad += (offset - dir * offset.dot(dir)).norm();
          n++;
        }
        sections[sec].radius = std::min(std::min(rad / n, max_radius), sections[par].radius);
      }
    }

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections[sec].ends.size() > 0)
    {
      BranchSection new_node;
      new_node.parent = (int)sec;
      new_node.roots = sections[sec].ends;
      new_node.max_distance_to_end = 0.0;
      for (auto &root : new_node.roots)
      {
        new_node.max_distance_to_end = std::max(new_node.max_distance_to_end, points[root].distance_to_end);
      }
      max_radius = radFromLength(new_node.max_distance_to_end, params);
      // constrain each new node's radius to be no larger than its parent radius. This is a reasonable constraint.
      new_node.radius = std::min(sections[sec].radius, max_radius);
      if (new_node.radius > 0.5 * params.min_diameter)  // if it is the first node, then we need a second node
      {
        sections[sec].children.push_back((int)sections.size());

        new_node.tip.setZero();
        for (auto &end : new_node.roots) 
        {
          new_node.tip += points[end].pos;
        }
        if (new_node.roots.size() > 0)
        {
          new_node.tip /= (double)new_node.roots.size();
        }
        sections.push_back(new_node);
      }
    }
  } // end of loop. We now have created all of the BranchSections
  
  // Now calculate the section ids for all of the points, for the segmented cloud
  for (size_t sec = 0; sec < sections.size(); sec++)
  {
    std::vector<int> nodes;
    if (sections[sec].ends.size() > 0)
    {
      for (auto &end : sections[sec].ends)
      {
        int node = points[end].parent;
        while (node != -1)
        {
          if (std::find(nodes.begin(), nodes.end(), node) != nodes.end())
          {
            break;
          }
          nodes.push_back(node);
          if (std::find(sections[sec].roots.begin(), sections[sec].roots.end(), node) != sections[sec].roots.end())
          {
            break;
          }
          node = points[node].parent;
        }
      }
    }
    std::vector<int> ends = sections[sec].ends;
    if (sections[sec].children.size() == 0)  // then also add any offshoots...
    {
      for (size_t i = 0; i < ends.size(); i++)
      {
        for (auto &c : children[ends[i]]) 
        {
          ends.push_back(c);
        }
      }
    }
    for (auto &node : nodes) 
    {
      section_ids[node] = (int)sec;
    }
    for (auto &end : ends) 
    {
      section_ids[end] = (int)sec;
    }
  }

  // for points without a section id (e.g. end points) 
  // look through their parents
  for (size_t i = 0; i < points.size(); i++)
  {
    if (section_ids[i] == -1)  // an unfound point
    {
      int j = points[i].parent;
      while (j != -1 && section_ids[j] == -1) 
      {
        j = points[j].parent;
      }
      if (j != -1)
      {
        section_ids[i] = section_ids[j];
      }
    }
  }
  // debug draw to rviz the set of cylinders
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<double> radii;
    for (auto &tree_node : sections)
    {
      if (tree_node.tip == Eigen::Vector3d(0, 0, 0))
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

  // generate local parent links. These are parent indices that are
  // independent per tree
  int num = 0;
  for (auto &section : sections)
  {
    if (section.parent >= 0 || section.children.empty())
      continue;
    num++;
    int child_id = 0;
    section.id = child_id++;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      sections[children[c]].id = child_id++;
      for (auto i : sections[children[c]].children) children.push_back(i);
    }
  }
  std::cout << num << " trees saved" << std::endl;

  Eigen::Vector3d min_bound(0, 0, 0), max_bound(0, 0, 0);
  // remove all sections with a root out of bounds, if we have gridded the cloud with an overlap
  if (params.grid_width)
  {
    double width = params.grid_width;
    cloud.calcBounds(&min_bound, &max_bound);
    Eigen::Vector3d mid = (min_bound + max_bound) / 2.0;
    Eigen::Vector2i inds(std::round(mid[0] / width), std::round(mid[1] / width));
    min_bound[0] = width * ((double)inds[0] - 0.5);
    min_bound[1] = width * ((double)inds[1] - 0.5);
    max_bound[0] = width * ((double)inds[0] + 0.5);
    max_bound[1] = width * ((double)inds[1] + 0.5);
    std::cout << "min bound: " << min_bound.transpose() << ", max bound: " << max_bound.transpose() << std::endl;

    // disable trees out of bounds
    for (auto &section : sections)
    {
      if (section.parent >= 0 || section.children.empty())
        continue;
      Eigen::Vector3d pos = section.tip;
      if (pos[0] < min_bound[0] || pos[0] > max_bound[0] || pos[1] < min_bound[1] || pos[1] > max_bound[1])
      {
        section.children.clear();  // make it a non-tree
      }
    }
  }

  // now colour the ray cloud based on the segmentation
  int j = -1;
  std::vector<int> root_segs(cloud.ends.size(), -1);
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    RGBA &colour = cloud.colours[i];
    if (cloud.rayBounded(i))
    {
      j++;
      int seg = section_ids[j];
      int root_id = points[j].root;
      root_segs[i] = root_id == -1 ? -1 : section_ids[root_id];

      if (!params.segment_branches)
      {
        if (root_id == -1)
        {
          colour.red = colour.green = colour.blue = 0;
          continue;
        }
        seg = section_ids[root_id];
      }
      if (seg == -1)
      {
        colour.red = colour.green = colour.blue = 0;
        continue;
      }
      convertIntToColour(seg, colour);
    }
    else
    {
      colour.red = colour.green = colour.blue = 0;
    }
  }

  if (params.grid_width)  // also remove edges from the segmented cloud
  {
    for (int i = (int)cloud.ends.size() - 1; i >= 0; i--)
    {
      if (!cloud.rayBounded(i))
        continue;
      Eigen::Vector3d pos = root_segs[i] == -1 ? cloud.ends[i] : sections[root_segs[i]].tip;

      if (pos[0] < min_bound[0] || pos[0] > max_bound[0] || pos[1] < min_bound[1] ||
          pos[1] > max_bound[1])  // nope, can't do this here!
      {
        cloud.starts[i] = cloud.starts.back();
        cloud.starts.pop_back();
        cloud.ends[i] = cloud.ends.back();
        cloud.ends.pop_back();
        cloud.colours[i] = cloud.colours.back();
        cloud.colours.pop_back();
        cloud.times[i] = cloud.times.back();
        cloud.times.pop_back();
      }
    }
  }
}

// save the structure to a text file
bool Trees::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << "# tree structure file:" << std::endl;
  ofs << "x,y,z,radius,parent_id,section_id" << std::endl; // simple format
  for (size_t sec = 0; sec < sections.size(); sec++)
  {
    const auto &section = sections[sec];
    if (section.parent >= 0 || section.children.empty()) // not a root section, so move on
    {
      continue;
    }
    ofs << section.tip[0] << "," << section.tip[1] << "," << section.tip[2] << "," << section.radius << ",-1," << sec;

    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      BranchSection &node = sections[children[c]];
      ofs << ", " << node.tip[0] << "," << node.tip[1] << "," << node.tip[2] << "," << node.radius << ","
          << sections[node.parent].id;
      ofs << ", " << children[c];
      for (auto i : sections[children[c]].children) 
      {
        children.push_back(i);
      }
    }
    ofs << std::endl;
  }
  return true;
}

}  // namespace ray
