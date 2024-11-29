// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include <nabo/nabo.h>
#include "rayclusters.h"

namespace ray
{
TreesParams::TreesParams()
  : max_diameter(0.9)
  , crop_length(1.0)
  , distance_limit(1.0)
  , height_min(2.0)
  , girth_height_ratio(0.12)
  , cylinder_length_to_width(4.0)
  , gap_ratio(0.016)
  , span_ratio(4.5)
  , gravity_factor(0.3)
  , grid_width(0.0)
  , segment_branches(false)
  , global_taper(0.012)
  , global_taper_factor(0.3)
{}

/// The main reconstruction algorithm
/// It is based on finding the shortest paths using Djikstra's algorithm, followed
/// by an agglomeration of paths, with repeated splitting from root to tips
Trees::Trees(Cloud &cloud, const Eigen::Vector3d &offset, const Mesh &mesh, const TreesParams &params, bool verbose)
{
  // firstly, get the full set of shortest paths from ground to tips, and the set of roots
  params_ = &params;

  std::vector<std::vector<int>> roots_list = getRootsAndSegment(
    points_, cloud, mesh, params_->max_diameter, params_->distance_limit, params_->height_min, params_->gravity_factor);

  // Now we want to convert these paths into a set of branch sections, from root to tips
  // splitting as we go up...

  // calculate the maximum distance to tip for every point in the cloud
  calculatePointDistancesToEnd();

  // then generate all the root sections
  generateRootSections(roots_list);

  // the set of child points will be useful later, so generate them now
  std::vector<std::vector<int>> children(points_.size());
  for (size_t i = 0; i < points_.size(); i++)
  {
    if (points_[i].parent != -1)
    {
      children[points_[i].parent].push_back(static_cast<int>(i));
    }
  }

  // first do a special case for all the trunks. This is where we estimate mean taper
  ray::Cloud debug_cloud;
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    double best_accuracy = -1.0;
    std::vector<int> nodes;  // all the points in the section

    // find height to tip by following the children:
    std::vector<int> list = sections_[sec_].roots;
    double max_height = -1e10;
    for (size_t j = 0; j < list.size(); j++)
    {
      max_height = std::max(max_height, points_[list[j]].pos[2]);
      list.insert(list.end(), children[list[j]].begin(), children[list[j]].end());
    }
    Eigen::Vector3d base = getRootPosition();
    double tree_height = std::max(0.01, max_height - base[2]);
    sections_[sec_].tree_height = tree_height;

    // Use a usr defined taper to control the height up the trunk to calculate the radius at
    double girth_height = params_->girth_height_ratio * tree_height; // sections_[sec_].max_distance_to_end;
    double estimated_radius = 1e10;
    double best_dist = 0.0;
    Eigen::Vector3d best_tip;
    
    std::vector<int> best_nodes;
    std::vector<int> best_ends;
    for (int j = 1; j<=3; j++)
    {
      double max_dist = girth_height * (double)j / 2.0; // range from 0.5 to 1.5 times the specified height
      nodes.clear();
      sections_[sec_].ends.clear();
      extractNodesAndEndsFromRoots(nodes, base, children, max_dist * 2.0/3.0, max_dist);
      if (nodes.size() < 2)
      {
        continue;
      }
      sections_[sec_].tip = calculateTipFromVertices(nodes);
      if (verbose)
      {
        for (auto &node: nodes)
        {
          debug_cloud.addRay(Eigen::Vector3d(0,0,0), points_[node].pos, 0.0, ray::RGBA(j==1 ? 255 : 0, j==2 ? 255:0, j==3 ? 255:0, 255));
        }
      }         
      if (removeDistantPoints(nodes))
      {
        sections_[sec_].tip = calculateTipFromVertices(nodes);
      }
      if (verbose)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip + Eigen::Vector3d(0,0,0.03), 0.0, ray::RGBA(255,255,0, 255));
      }
      // shift to cylinder's centre
      Eigen::Vector3d up(0,0,1);
      sections_[sec_].tip += vectorToCylinderCentre(nodes, up);
      // now find the segment radius
      double accuracy;
      double radius = estimateCylinderRadius(nodes, up, accuracy);

      ray::RGBA col(j==1 ? 255 : 127, j==2 ? 255:127, j==3 ? 255:127, 255);
      if (verbose)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip, 0.0, col);
        for (double ang = 0; ang < 2.0*ray::kPi; ang += 0.1)
        {
          debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip + radius * Eigen::Vector3d(std::sin(ang), std::cos(ang),0), 0.0, col);
        }
      }
      if (radius < estimated_radius)
      {
        best_accuracy = accuracy;
        estimated_radius = radius;
        best_dist = max_dist;
        best_tip = sections_[sec_].tip;
        best_nodes = nodes;
        best_ends = sections_[sec_].ends;
      }
    }
    if (best_dist == 0.0)
    {
      std::cout << "warning: could not find any points on trunk " << sec_ << " at " << base.transpose() << " so removing the whole section" << std::endl;
      sections_[sec_].tip = base + Eigen::Vector3d(0,0,0.01);
      sections_[sec_].total_weight = 1e-10;
      sections_[sec_].ends.clear(); // so this trunk is not ever used
      continue; 
    }    
    sections_[sec_].tip = best_tip;
    sections_[sec_].ends = best_ends;
    nodes = best_nodes;
    if (sections_[sec_].split_count < 2)
    {
      double thickness = best_dist; 
      bool points_removed = false;
      double gap = params_->gap_ratio * sections_[sec_].max_distance_to_end; // gap threshold for splitting
      double span = params_->span_ratio * estimated_radius; // span threshold for splitting
      std::vector<std::vector<int>> clusters = findPointClusters(base, points_removed, thickness, span, gap);

      if (clusters.size() > 1 || (points_removed && clusters.size() > 0))  // a bifurcation (or an alteration)
      {
        sections_[sec_].split_count++;
        bifurcate(clusters, thickness, children, true, true);
        sec_--;
        continue;
      }
    }
    if (verbose)
    {
      for (auto &node: sections_[sec_].ends)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), points_[node].pos + Eigen::Vector3d(0,0,0.02), 0.0, ray::RGBA(255, 0, 255, 255));
      }
      for (double ang = 0; ang < 2.0*ray::kPi; ang += 0.1)
      {
        uint8_t shade = 255; // (uint8_t)(best_accuracy * 255.0);
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), best_tip + (estimated_radius + 0.01) * Eigen::Vector3d(std::sin(ang), std::cos(ang),0), 0.0, ray::RGBA(shade,shade,shade,255));
      }
    }

    nodes.clear();
    sections_[sec_].ends.clear();
    extractNodesAndEndsFromRoots(nodes, base, children, 0.0, best_dist/2.0); // make it lower
    estimateCylinderTaper(estimated_radius, best_accuracy, false); // update the expected taper
  }
  if (verbose)
  {
    debug_cloud.translate(offset);
    debug_cloud.save("debug.ply");
  }

  // now trace from root tree nodes upwards, getting node centroids
  // create new BranchSections as we go
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    const int par = sections_[sec_].parent;
    if (par == -1)
    {
      // now add the single child for this particular tree node, assuming there are still ends
      if (sections_[sec_].ends.size() > 0)
      {
        addChildSection();
        for (const auto &i: sections_[sec_].roots)
        {
          sections_[sec_].tip[2] = std::min(sections_[sec_].tip[2], points_[i].pos[2]);
        }
      }      
      else
      {
        std::cout << "weird, a trunk without end points! " << sec_ << std::endl;
      }
      continue;
    }
    if (!(sec_ % 10000))
    {
      std::cout << "generating segment " << sec_ << std::endl;
    }

    // this branch can come to an end as it is now too small
    if (sections_[sec_].max_distance_to_end < params_->crop_length)
    {
      setBranchTip();
      continue;
    }

    std::vector<int> nodes;  // all the points in the section
    bool extract_from_ends = sections_[sec_].ends.size() > 0;
    // if the branch section has no end points recorded, then we need to examine this branch to
    // find end points and potentially branch (bifurcate)
    if (!extract_from_ends)
    {
      Eigen::Vector3d base = getRootPosition();
      double span_rad = radius(sections_[sec_]); 
      double thickness = params_->cylinder_length_to_width * span_rad;
      extractNodesAndEndsFromRoots(nodes, base, children, 0.0, thickness);
      
      bool points_removed = false;
      double gap = params_->gap_ratio * sections_[sec_].max_distance_to_end; // gap threshold for splitting
      double span = params_->span_ratio * span_rad; // span thershold for splitting
      std::vector<std::vector<int>> clusters = findPointClusters(base, points_removed, thickness, span, gap);

      if (clusters.size() > 1 || (points_removed && clusters.size() > 0))  // a bifurcation (or an alteration)
      {
        extract_from_ends = true; // don't trust the found nodes as it is now two separate tree nodes
        bool add_offshoots = true; // par != -1 && sections_[par].parent != -1;  // when this is false it can lead to whole branches missing.
        // if points have been removed then this only resets the current section's points
        // otherwise it creates new branch sections_ for each cluster and adds to the end of the sections_ list
        bifurcate(clusters, thickness, children, false, add_offshoots);
      }
    }

    if (extract_from_ends) // we have split the ends, so we need to extract the set of nodes in a backwards manner
    {
      extractNodesFromEnds(nodes);
    }
    // estimate the section's tip (the centre of the cylinder of points)
    sections_[sec_].tip = calculateTipFromVertices(nodes);
    // get section's direction
    Eigen::Vector3d dir = par >= 0 ? (sections_[sec_].tip - sections_[par].tip).normalized() : Eigen::Vector3d(0, 0, 1);
    // shift to cylinder's centre
    sections_[sec_].tip += vectorToCylinderCentre(nodes, dir);
    // re-estimate direction
    dir = par >= 0 ? (sections_[sec_].tip - sections_[par].tip).normalized() : Eigen::Vector3d(0, 0, 1);
    // now find the segment radius
    double accuracy = 0.0;
    double rad = estimateCylinderRadius(nodes, dir, accuracy);
    // and estimate taper
    estimateCylinderTaper(rad / sections_[sec_].radius_scale, accuracy, extract_from_ends);

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections_[sec_].ends.size() > 0)
    {
      addChildSection();
    }
  }  // end of loop. We now have created all of the BranchSections

  // Now calculate the section ids for all of the points, for the segmented cloud
  std::vector<int> section_ids(points_.size(), -1);
  calculateSectionIds(section_ids, children);

  generateLocalSectionIds();

  Eigen::Vector3d min_bound(0, 0, 0), max_bound(0, 0, 0);
  // remove all sections with a root out of bounds, if we have gridded the cloud with an overlap
  if (params_->grid_width)
  {
    removeOutOfBoundSections(cloud, min_bound, max_bound, offset);
  }

  std::vector<int> root_segs(cloud.ends.size(), -1);
  // now colour the ray cloud based on the segmentation
  segmentCloud(cloud, root_segs, section_ids);

  if (params_->grid_width)  // also remove rays from the segmented cloud
  {
    removeOutOfBoundRays(cloud, min_bound, max_bound, root_segs);
  }

  std::cout << "cloud's estimated mean taper ratio (diameter / length): " << 2.0 * forest_taper_ / forest_weight_ << std::endl;
}

// If 1 tree plus some low lying foliage is in the node, then it won't be split if the foliage doesn't reach to the 
// top of the section (so doesn't have an end point), but it can create a very large radius estimate. 
// here we remove points that are too far from any end position. 
// this allows us to use most nodes to get an accurate radius estimate, but still uses only the end nodes to decide whether
// to split
bool Trees::removeDistantPoints(std::vector<int> &nodes)
{
  Eigen::Vector3d up(0,0,1);
  double max_rad = params_->gap_ratio * sections_[sec_].max_distance_to_end; // gap threshold for splitting
  double max_rad_sqr = max_rad * max_rad;
  bool found = false;
  for (int k = 0; k<(int)nodes.size(); k++)
  {
    if (nodes.size() <= 6)
    {
      break;
    }
    bool all_distant = true;
    for (auto &end: sections_[sec_].ends)
    {
      Eigen::Vector3d dif = points_[nodes[k]].pos - points_[end].pos;
      dif[2] = 0.0;
      if (dif.squaredNorm() < max_rad_sqr)
      {
        all_distant = false;
        break;
      }
    }
    if (all_distant)
    {
      nodes[k] = nodes.back();
      nodes.pop_back();
      k--;
      found = true;
    }
  }
  return found;
}

double Trees::meanTaper(const BranchSection &section) const
{
  int root = section.root;
  double mean_taper = params_->global_taper ? params_->global_taper : (forest_taper_ / forest_weight_);
  double mean_weight = forest_weight_squared_ / forest_weight_;
  double taper = sections_[root].total_taper / sections_[root].total_weight;
  double weight = sections_[root].total_weight;
  double blend = params_->global_taper_factor * params_->global_taper_factor * params_->global_taper_factor;
  mean_weight *= blend;
  weight *= 1.0 - blend;

  double result = (mean_taper * mean_weight + taper*weight) / (mean_weight + weight);
  if (!(result == result))
  {
    std::cout << "bad taper estimate" << std::endl;
    std::cout << "root: " << root << " estimated taper: " << result << ", mt: " << mean_taper << ", mw: " << mean_weight << ", t: " << taper << ", w: " << weight << ", blend: " << blend << std::endl;
  }
  return result;
}

double Trees::radius(const BranchSection &section) const 
{ 
  return sections_[section.root].tree_height * meanTaper(section) * section.radius_scale; // scaled down from root's estimated radius
}

void Trees::calculatePointDistancesToEnd()
{
  for (size_t i = 0; i < points_.size(); i++)
  {
    int id = static_cast<int>(i);
    int parent = points_[i].parent;
    points_[i].visited = false;  // reset the visited flag while we're here
    while (parent != -1)
    {
      const double end_dist =
        points_[id].distance_to_end + (points_[id].distance_to_ground - points_[parent].distance_to_ground);
      if (points_[parent].distance_to_end > end_dist)
      {
        break;
      }
      points_[parent].distance_to_end = end_dist;
      id = parent;
      parent = points_[id].parent;
    }
  }
}

// generate the initial branch sections, which are the root locations, defined by @c roots_list
// we also fill in the @c section_ids, which record which section each point belongs to
void Trees::generateRootSections(const std::vector<std::vector<int>> &roots_list)
{
  // first we initialise to one branch section per root
  for (size_t i = 0; i < roots_list.size(); i++)
  {
    BranchSection base;
    base.roots = roots_list[i];
    base.root = static_cast<int>(sections_.size());
    base.radius_scale = 1.0;
    sections_.push_back(base);
  }

  // we can now find the maximum distance to tip per root section, which allows us to estimate the
  // section radius for all the root sections.
  for (int i = static_cast<int>(sections_.size()) - 1; i >= 0; i--)
  {
    for (auto &root : sections_[i].roots)
    {
      sections_[i].max_distance_to_end = std::max(sections_[i].max_distance_to_end, points_[root].distance_to_end);
    }
  }
}

void Trees::setBranchTip()
{
  // before we finish, we need to calculate the tip of this branch, from its end points
  sections_[sec_].tip.setZero();
  for (auto &i : sections_[sec_].ends)
  {
    sections_[sec_].tip += points_[i].pos;
  }
  size_t sum = sections_[sec_].ends.size();
  // if it doesnt have end points then we'll have to use the root points
  if (sections_[sec_].ends.empty())
  {
    for (auto &i : sections_[sec_].roots)
    {
      sections_[sec_].tip += points_[i].pos;
    }
    sum = sections_[sec_].roots.size();
  }
  if (sum > 0)
  {
    sections_[sec_].tip /= static_cast<double>(sum);
  }
}

Eigen::Vector3d Trees::getRootPosition() const
{
  // get a base position for the section
  Eigen::Vector3d base(0, 0, 0);
  for (auto &root : sections_[sec_].roots)
  {
    base += points_[root].pos;
  }
  base /= static_cast<double>(sections_[sec_].roots.size());
  return base;
}


// find the node and end points for this section, from the root points
// note that the nodes contain everything between min_dist and max_dist and also the adjacent nodes either side
// i.e. a root node and an end node. This means the minimum number of nodes is 2
void Trees::extractNodesAndEndsFromRoots(std::vector<int> &nodes, const Eigen::Vector3d &base,
                                         const std::vector<std::vector<int>> &children, double min_dist, double max_dist)
{
  const int par = sections_[sec_].parent;
  std::vector<int> iteration_nodes = sections_[sec_].roots;

  // 2. find all the points (and the end points) for this section:
  for (unsigned int index = 0; index < iteration_nodes.size(); index++)
  {
    // we cannot use the for (int i: nodes) syntax above, as we are pushing into nodes as we go
    const int i = iteration_nodes[index];
    bool use_parent = false;
    for (auto &child : children[i])
    {
      const double dist = par == -1 ? points_[child].pos[2] - base[2] : (points_[child].pos - base).norm();
      if (dist > min_dist)
      {
        if (index < sections_[sec_].roots.size())
        {
          if (!use_parent)
          {
            nodes.push_back(i);
          }
          use_parent = true;
        }
        nodes.push_back(child);
      }
      if (dist < max_dist)  // in same slot, so accumulate
      {
        iteration_nodes.push_back(child);  // so we recurse on this child too
      }
      else
      {
        sections_[sec_].ends.push_back(child);
      }
    }
  }
}

// find clusters of points from the root points up the shortest paths, up to the cylinder length
std::vector<std::vector<int>> Trees::findPointClusters(const Eigen::Vector3d &base, bool &points_removed,
                              double thickness, double span, double gap)
{
  const int par = sections_[sec_].parent;
  std::vector<int> &all_ends = sections_[sec_].ends;

  // convert to a structure that is better for the cluster function
  std::vector<Eigen::Vector3d> ps;
  std::vector<int> v_indices;
  ps.reserve(all_ends.size());
  v_indices.reserve(all_ends.size());

  // cluster end points to find if we have separate branches
  // first, get interpolated edge points_. i.e. interpolate between two connected points inside and outside
  // the branch section's cylinder, so that the set edge_pos are right on the top boundary of the section
  for (auto &j : all_ends)
  {
    const double dist1j = par == -1 ? points_[j].pos[2] - base[2] : (points_[j].pos - base).norm();
    const double dist0j =
      par == -1 ? points_[points_[j].parent].pos[2] - base[2] : (points_[points_[j].parent].pos - base).norm();
    const double blendj = dist1j == dist0j ? 0.0 : (thickness - dist0j) / (dist1j - dist0j);
    Eigen::Vector3d edge_pos = points_[points_[j].parent].pos * (1.0 - blendj) + points_[j].pos * blendj;
    ps.push_back(edge_pos);
    v_indices.push_back(j);
  }

  // cluster these end points based on two separation criteria (gap_ratio and span_ratio)
  std::vector<std::vector<int>> clusters;
  generateClusters(clusters, ps, gap, span);
  // adjust back to global ids
  for (auto &cluster : clusters)
  {
    for (auto &id : cluster)
    {
      id = v_indices[id];
    }
  }

  points_removed = false;
  if (par == -1)  // if this is the root section
  {
    // then remove children that are smaller than the minimum tree height
    for (int i = static_cast<int>(clusters.size()) - 1; i >= 0; i--)
    {
      double max_dist = 0.0;
      for (auto &end : clusters[i])
      {
        max_dist = std::max(max_dist, points_[end].distance_to_end);
      }
      if (max_dist < params_->height_min)
      {
        clusters[i] = clusters.back();
        clusters.pop_back();
        points_removed = true;
      }
    }
  }
  return clusters;
}

// split into multiple branches and add as new branch sections to the end of the sections_ list
// that is being iterated through.
void Trees::bifurcate(const std::vector<std::vector<int>> &clusters, double thickness, std::vector<std::vector<int>> &children, bool clip_tree, bool add_offshoots)
{
  const int par = sections_[sec_].parent;
  // find the maximum distance (to tip) for each cluster
  std::vector<double> max_distances(clusters.size(), 0.0);
  double maxmax = -1;
  int maxi = -1;
  for (size_t i = 0; i < clusters.size(); i++)
  {
    for (auto &end : clusters[i])
    {
      max_distances[i] = std::max(max_distances[i], points_[end].distance_to_end);
    }
    max_distances[i] += thickness;
    if (max_distances[i] > maxmax)
    {
      maxmax = max_distances[i];
      maxi = static_cast<int>(i);
    }
  }

  double total_area = 0.0;
  for (auto &dist: max_distances)
  {
    total_area += dist * dist;
  }

  // if clip_tree then we need to change the root nodes to be only those ones that
  // each path goes to... 
  
  // set the current branch section to the cluster with the
  // maximum distance to tip (i.e. the longest branch).
  sections_[sec_].ends = clusters[maxi];
  std::vector<int> all_nodes;
  std::vector<int> main_roots;

  if (clip_tree)
  {
    for (auto &end : sections_[sec_].ends)
    {
      int node = points_[end].parent;
      while (node != -1)
      {
        if (std::find(all_nodes.begin(), all_nodes.end(), node) != all_nodes.end())
        {
          break;
        }
        all_nodes.push_back(node);  // fill in the nodes in this branch section
        if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
        {
          main_roots.push_back(node);
          break;
        }
        node = points_[node].parent;
      }
    }
  }

  sections_[sec_].max_distance_to_end = max_distances[maxi];

  if (add_offshoots)
  {
    // for all other clusters, add new sections to the list...
    for (size_t i = 0; i < clusters.size(); i++)
    {
      if (static_cast<int>(i) == maxi)
      {
        continue;
      }
      BranchSection new_node = sections_[sec_];
      new_node.max_distance_to_end = max_distances[i];
      new_node.radius_scale *= max_distances[i] / std::sqrt(total_area);
      if (par == -1)
      {
        new_node.radius_scale = 1.0;
      }
      if (new_node.max_distance_to_end > params_->crop_length)  // but only add if they are large enough
      {
        // we only specify the end points at this stage. They will therefore enter the
        // extract_from_ends block below when the sec_ for loop gets to their section
        new_node.ends = clusters[i];

        if (clip_tree)
        {
          std::vector<int> my_nodes;
          std::vector<int> my_roots;
          for (auto &end : new_node.ends)
          {
            int node = points_[end].parent;
            int child = end;
            while (node != -1)
            {
              if (std::find(my_nodes.begin(), my_nodes.end(), node) != my_nodes.end()) 
              {
                break; // just hit my own tree
              }
              if (std::find(all_nodes.begin(), all_nodes.end(), node) != all_nodes.end())
              {
                // hit another tree, so make a root here
                my_roots.push_back(child);
                // more than this, we need to unlink the children bit.
                auto &list = children[node];

                int find_count = 0;
                for (int j = 0; j<(int)list.size(); j++)
                {
                  if (list[j] == child)
                  {
                    find_count++;
                    list[j] = list.back();
                    list.pop_back();
                    j--;
                  }
                }
                points_[child].parent = -1;

                // we need to tell all of the children what the new root is
                std::vector<int> child_list = {child};
                for (size_t j = 0; j<child_list.size(); j++)
                {
                  points_[child_list[j]].root = child;
                  child_list.insert(child_list.end(), children[child_list[j]].begin(), children[child_list[j]].end());
                }
                break;
              }
              my_nodes.push_back(node);  // fill in the nodes in this branch section
              if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
              {
                my_roots.push_back(node);
                break;
              }
              child = node;
              node = points_[node].parent;
            }
          }
          all_nodes.insert(all_nodes.end(), my_nodes.begin(), my_nodes.end());
          new_node.roots = my_roots;
        }

        if (par != -1)
        {
          sections_[par].children.push_back(static_cast<int>(sections_.size()));
        }
        if (par == -1)
        {
          new_node.root = (int)sections_.size();
        }
        sections_.push_back(new_node);
      }
    }
  }
  if (clip_tree)
  {
    sections_[sec_].roots = main_roots;
  }
  if (par > -1)
  {
    sections_[sec_].radius_scale *= max_distances[maxi] / std::sqrt(total_area); // reduce branch radius due to the split
  }
}

// find the nodes between the section end points and the section root points
void Trees::extractNodesFromEnds(std::vector<int> &nodes)
{
  nodes = sections_[sec_].ends;
  for (auto &end : sections_[sec_].ends)
  {
    int node = points_[end].parent;
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
      nodes.push_back(node);  // fill in the nodes in this branch section
      if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
      {
        break;
      }
      node = points_[node].parent;
    }
  }
}

// use the list of Vertex IDs to estimate a tip location
Eigen::Vector3d Trees::calculateTipFromVertices(const std::vector<int> &vertices) const
{
  Eigen::Vector3d tip(0, 0, 0);
  if (vertices.empty())
  {
    std::cerr << "error: there shouldn't be empty nodes at this point in the processing" << std::endl;
  }
  for (auto &i : vertices)
  {
    tip += points_[i].pos;
  }
  if (!vertices.empty())
  {
    tip /= static_cast<double>(vertices.size());
  }

  return tip;
}

Eigen::Vector3d Trees::vectorToCylinderCentre(const std::vector<int> &nodes, const Eigen::Vector3d &dir) const
{
#define REAL_CENTROID  // finds a new centroid that is robust to branches scanned from a single side
#if defined REAL_CENTROID
  Eigen::Vector3d mean_p(0, 0, 0);
  std::vector<Eigen::Vector3d> ps;
  const Eigen::Vector3d vec(1, 2, 3);
  // obtain two orthogonal planes to the section's direction vector
  const Eigen::Vector3d ax1 = dir.cross(vec).normalized();
  const Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
  double n = 0;
  // project points into a local space parabola
  for (auto &i : nodes)
  {
    const Eigen::Vector3d pos = points_[i].pos - sections_[sec_].tip;
    const Eigen::Vector2d offset(ax1.dot(pos), ax2.dot(pos));
    const Eigen::Vector2d xy = offset;
    const double l2 = xy.squaredNorm();
    const Eigen::Vector3d point(xy[0], xy[1], 0.5 * l2);  // a paraboloid that has gradient 1 at 1
    ps.push_back(point);
    mean_p += point;
    n++;
  }
  mean_p /= n;
  double approx_rad_sqr = 2.0 * mean_p[2];
  if (n > 5)  // assuming there are sufficient points for a resonable guess
  {
    // accumulation structure for plane least squares fitting
    struct Acc
    {
      Acc() { x2 = y2 = xy = xz = yz = 0; }
      double x2, y2, xy, xz, yz;
    };
    Acc plane;
    for (const auto &p : ps)
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
      const double A = (plane.xz * plane.y2 - plane.yz * plane.xy) / (plane.x2 * plane.y2 - plane.xy * plane.xy);
      const double B = (plane.yz - A * plane.xy) / plane.y2;

      Eigen::Vector2d shift(A, B);
      const double shift2 = shift.squaredNorm();
      if (shift2 > approx_rad_sqr)  // don't shift more than one radius each iteration, for safety
      {
        shift *= std::sqrt(approx_rad_sqr / shift2);
      }

      // apply the plane parameters as a world-space shift in the branch section tip position
      return ax1 * shift[0] + ax2 * shift[1];
    }
  }
  else if (n > 0)  // if only a few points then use the mean instead
  {
    return ax1 * mean_p[0] + ax2 * mean_p[1];
  }
#endif
  return Eigen::Vector3d(0, 0, 0);
}

double Trees::estimateCylinderRadius(const std::vector<int> &nodes, const Eigen::Vector3d &dir, double &accuracy)
{
  double rad = 0.0;
  // get the mean radius
  double power = 0.25; // when 1 this is usual radius, but lower powers reduce outliers due to foliage
  std::vector<double> norms;
  norms.reserve(nodes.size());
  if (params_->use_rays)
  {
    for (auto &node : nodes)
    {
      Eigen::Vector3d offset = points_[node].pos - sections_[sec_].tip;
      offset -= dir * offset.dot(dir); // flatten    
      Eigen::Vector3d ray = points_[node].start - points_[node].pos;
      ray -= dir * ray.dot(dir); // flatten
      offset += ray*std::max(0.0, std::min(-offset.dot(ray)/ray.dot(ray), 1.0)); // move to closest point
      norms.push_back(offset.norm());
    }
  }
  else
  {
    for (auto &node : nodes)
    {
      Eigen::Vector3d offset = points_[node].pos - sections_[sec_].tip;
      offset -= dir * offset.dot(dir); // flatten    
      norms.push_back(offset.norm());
    }
  }

  for (auto &norm: norms)
  {
    rad += std::pow(norm, power);
  }
  const double eps = 1e-5; // prevent division by 0
  rad /= (double)nodes.size() + eps;
  rad = std::pow(rad, 1.0/power);
  double e = 0.0;
  for (auto &norm : norms)
  {
    e += std::abs(norm - rad);
  }
  e /= (double)nodes.size() + eps;
#define NEW_ACCURACY
#if defined NEW_ACCURACY
  const double sensor_noise = 0.02; // this prevents cylinders with a small number of very accurate points from totally dominating
  accuracy = rad / (e + sensor_noise);
#else // this one goes down to 0, which could be bad if all cylinder points were low accuracy
  accuracy = std::max(eps, rad - 2.0*e) / std::max(eps, rad); // so even distribution is accuracy 0, and cylinder shell is accuracy 1 
  accuracy *= std::min((double)nodes.size() / 3.0, 1.0);
#endif
  accuracy = std::max(accuracy, eps);
  rad = std::max(rad, eps);

  return rad;
}

void Trees::estimateCylinderTaper(double radius, double accuracy, bool extract_from_ends)
{
  int par = sections_[sec_].parent;
  int root = sections_[sec_].root;

  double l = sections_[sec_].radius_scale * sections_[root].tree_height;
  double L = sections_[root].tree_height;

  double junction_weight = 1.0;
  if (par > -1)
  {
    junction_weight = extract_from_ends ? 0.25 : 1.0; // extracting from ends means we are less confident about the quality
  }

  double weight = l*l*l * accuracy * junction_weight;
 // weight *= weight; // preference the strongest weight sections
  double taper = (radius/L) * weight;

  forest_taper_ += taper;
  forest_weight_ += weight;
  forest_weight_squared_ += weight*weight;

  sections_[root].total_taper += taper;
  sections_[root].total_weight += weight;
  if (sections_[root].total_weight == 0.0)
  {
    std::cout << "bad section weight: " << l << "," << accuracy << "," << junction_weight << std::endl;
    std::cout << "sec: " << sec_ << ", root: " << root << std::endl;
  }

  sections_[sec_].taper = taper;
  sections_[sec_].weight = weight;
  sections_[sec_].len = l*l*l;
  sections_[sec_].accuracy = accuracy;
  sections_[sec_].junction_weight = junction_weight;
}

// add a child section to continue reconstructing the tree segments
void Trees::addChildSection()
{
  BranchSection new_node;
  new_node.parent = static_cast<int>(sec_);
  new_node.root = sections_[sec_].root;
  new_node.roots = sections_[sec_].ends;
  new_node.taper = sections_[sec_].taper;
  new_node.weight = sections_[sec_].weight;
  new_node.radius_scale = sections_[sec_].radius_scale;
  
  for (auto &root : new_node.roots)
  {
    new_node.max_distance_to_end = std::max(new_node.max_distance_to_end, points_[root].distance_to_end);
  }
  if (new_node.max_distance_to_end > params_->crop_length)  
  {
    sections_[sec_].children.push_back(static_cast<int>(sections_.size()));
    new_node.tip = calculateTipFromVertices(new_node.roots);
    sections_.push_back(new_node);
  }
}

// calculate what section every point belongs to
void Trees::calculateSectionIds(std::vector<int> &section_ids, const std::vector<std::vector<int>> &children)
{
  // first we initialise to one branch section per root
  for (size_t i = 0; i < sections_.size(); i++)
  {
    if (sections_[i].parent == -1)
    {
      for (auto &root : sections_[i].roots)
      {
        section_ids[root] = static_cast<int>(i);
      }
    }
  }
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    std::vector<int> nodes;
    for (auto &end : sections_[sec_].ends)
    {
      int node = points_[end].parent;
      while (node != -1)
      {
        if (std::find(nodes.begin(), nodes.end(), node) != nodes.end())
        {
          break;
        }
        nodes.push_back(node);
        if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) !=
            sections_[sec_].roots.end())
        {
          break;
        }
        node = points_[node].parent;
      }
    }
    std::vector<int> ends = sections_[sec_].ends;
    if (sections_[sec_].children.size() == 0)  // then also add any offshoots...
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
      section_ids[node] = static_cast<int>(sec_);
    }
    for (auto &end : ends)
    {
      section_ids[end] = static_cast<int>(sec_);
    }
  }

  // for points without a section id (e.g. end points)
  // look through their parents
  for (size_t i = 0; i < points_.size(); i++)
  {
    if (section_ids[i] == -1)  // an unfound point
    {
      int j = points_[i].parent;
      while (j != -1 && section_ids[j] == -1)
      {
        j = points_[j].parent;
      }
      if (j != -1)
      {
        section_ids[i] = section_ids[j];
      }
    }
  }
}

// generate local parent links. These are parent indices that are
// independent per tree
void Trees::generateLocalSectionIds()
{
  int num = 0;
  for (auto &section : sections_)
  {
    if (section.parent >= 0 || section.children.empty())
    {
      continue;
    }
    num++;
    int child_id = 0;
    section.id = child_id++;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      sections_[children[c]].id = child_id++;
      for (auto i : sections_[children[c]].children)
      {
        children.push_back(i);
      }
    }
  }
  std::cout << num << " trees saved" << std::endl;
}

void Trees::removeOutOfBoundSections(const Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound, const Eigen::Vector3d &offset)
{
  const double width = params_->grid_width;
  cloud.calcBounds(&min_bound, &max_bound);
  const Eigen::Vector3d mid = (min_bound + max_bound)/2.0 + offset;
  const Eigen::Vector2d inds(std::round(mid[0] / width), std::round(mid[1] / width));
  min_bound[0] = width * (inds[0] - 0.5) - offset[0];
  min_bound[1] = width * (inds[1] - 0.5) - offset[1];
  max_bound[0] = width * (inds[0] + 0.5) - offset[0];
  max_bound[1] = width * (inds[1] + 0.5) - offset[1];
  std::cout << "min bound: " << (min_bound+offset).transpose() << ", max bound: " << (max_bound+offset).transpose() << std::endl;

  // disable trees out of bounds
  for (auto &section : sections_)
  {
    if (section.parent >= 0 || section.children.empty())
    {
      continue;
    }
    const Eigen::Vector3d pos = section.tip;
    if (pos[0] < min_bound[0] || pos[0] > max_bound[0] || pos[1] < min_bound[1] || pos[1] > max_bound[1])
    {
      section.children.clear();  // make it a non-tree
    }
  }
}

// colour the cloud by tree id, or by branch segment id
void Trees::segmentCloud(Cloud &cloud, std::vector<int> &root_segs, const std::vector<int> &section_ids)
{
  contiguous_section_ids_.resize(sections_.size(), -1); // these are different to the root section IDs as they exclude empty trees
  int num_trees = 0;
  bool first_non_root = false;
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent < 0 && first_non_root)
    {
      std::cout << "Error: root sections after some non-roots. ID: " << sec << std::endl;
    }
    if (section.parent >= 0)
      first_non_root = true;
    else if (section.children.empty())  // not a root section, so move on
    {
      continue;
    }
    contiguous_section_ids_[sec] = num_trees++;
  }

  int j = -1;
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    RGBA &colour = cloud.colours[i];
    if (cloud.rayBounded(i))
    {
      j++;
      int seg = section_ids[j];
      const int root_id = points_[j].root;
      root_segs[i] = root_id == -1 ? -1 : section_ids[root_id];

      // This block places a cone of gradient 2 around the trunk direction, to remove the square of initial ground points
      if (seg != -1 && sections_[seg].parent == -1)
      {
        Eigen::Vector3d dir(0,0,1e-10);
        for (auto &child : sections_[seg].children)
        {
          dir += sections_[child].tip - sections_[seg].tip;
        }
        dir.normalize();
        const double grad = 2.0; // larger cuts out a steeper (narrower) cone
        Eigen::Vector3d base = sections_[seg].tip - grad*dir*radius(sections_[seg]);
        Eigen::Vector3d dif = cloud.ends[i] - base;
        double h = dif.dot(dir);
        double w = (dif - dir*h).norm();
        if (grad*w > h)
        {
          colour.red = colour.green = colour.blue = 0;
          continue;
        }
      }

      if (!params_->segment_branches)
      {
        if (root_id == -1)
        {
          colour.red = colour.green = colour.blue = 0;
          continue;
        }
        seg = root_segs[i];
      }
      if (seg == -1)
      {
        colour.red = colour.green = colour.blue = 0;
        continue;
      }
      convertIntToColour(contiguous_section_ids_[seg], colour);
    }
    else
    {
      colour.red = colour.green = colour.blue = 0;
    }
  }
}

// remove rays from the ray cloud where the end points are out of bounds
void Trees::removeOutOfBoundRays(Cloud &cloud, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound,
                                 const std::vector<int> &root_segs)
{
  for (int i = static_cast<int>(cloud.ends.size()) - 1; i >= 0; i--)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    const Eigen::Vector3d pos = root_segs[i] == -1 ? cloud.ends[i] : sections_[root_segs[i]].tip;

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

// save the structure to a text file
bool Trees::save(const std::string &filename, const Eigen::Vector3d &offset, bool verbose) const
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << std::setprecision(4) << std::fixed;
  ofs << "# Tree file. Optional per-tree attributes (e.g. 'height,crown_radius, ') followed by 'x,y,z,radius' and any additional per-segment attributes:" << std::endl;
  ofs << "x,y,z,radius,parent_id,section_id"; // simple format
  if (verbose)
  {
    ofs << ",weight,len,accuracy,junction_weight";
  }
  ofs << std::endl;  
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent >= 0 || section.children.empty())  // not a root section, so move on
    {
      continue;
    }
    ofs << section.tip[0]+offset[0] << "," << section.tip[1]+offset[1] << "," << section.tip[2]+offset[2] << "," << radius(section) << ",-1," << contiguous_section_ids_[sec];
    if (verbose)
    {
      ofs << "," << section.weight << "," << section.len << "," << section.accuracy << "," << section.junction_weight;
    }
    int root = section.root;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      const BranchSection &node = sections_[children[c]];
      if (node.root != root)
      {
        std::cout << "bad format: " << node.root << " != " << root << std::endl;
      }
      ofs << ", " << node.tip[0]+offset[0] << "," << node.tip[1]+offset[1] << "," << node.tip[2]+offset[2] << "," << radius(node) << "," << sections_[node.parent].id << "," << contiguous_section_ids_[children[c]];
      if (verbose)
      {
        ofs << "," << node.weight << "," << node.len << "," << node.accuracy << "," << node.junction_weight;
      }
      for (auto i : sections_[children[c]].children)
      {
        children.push_back(i);
      }
    }
    ofs << std::endl;
  }
  return true;
}

}  // namespace ray
