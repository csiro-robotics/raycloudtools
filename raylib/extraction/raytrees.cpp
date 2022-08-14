// Copyright (c) 2022
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
TreesParams::TreesParams()
  : max_diameter(0.9)
  , min_diameter(0.02)
  , distance_limit(1.0)
  , height_min(2.0)
  , length_to_radius(140.0)
  , cylinder_length_to_width(4.0)
  , gap_ratio(2.5)
  , span_ratio(4.5)
  , gravity_factor(0.3)
  , radius_exponent(0.67)
  , linear_range(3.0)
  , grid_width(0.0)
  , segment_branches(false)
{}

/// The main reconstruction algorithm
/// It is based on finding the shortest paths using Djikstra's algorithm, followed
/// by an agglomeration of paths, with repeated splitting from root to tips
Trees::Trees(Cloud &cloud, const Mesh &mesh, const TreesParams &params, bool verbose)
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

  // now trace from root tree nodes upwards, getting node centroids
  // create new BranchSections as we go
  for (sec_ = 0; sec_ < sections_.size(); sec_++)
  {
    if (!(sec_ % 10000))
    {
      std::cout << "generating segment " << sec_ << std::endl;
    }
    const int par = sections_[sec_].parent;

    // any estimate of branch radius is bounded from above by an estimate from the branch's length
    // estimates from length are more reliable on small branches, and less reliable on large, well
    // observed trunks.
    max_radius_ = radFromLength(sections_[sec_].max_distance_to_end);

    // this branch can come to an end as it is now too small
    if (max_radius_ < 0.5 * params_->min_diameter)
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
      extractNodesAndEndsFromRoots(nodes, base, children);
      bool points_removed = false;
      std::vector<std::vector<int>> clusters = findPointClusters(base, points_removed);

      if (clusters.size() > 1 || points_removed)  // a bifurcation (or an alteration)
      {
        extract_from_ends = true;
        nodes.clear();  // don't trust the found nodes as it is now two separate tree nodes
        // if points have been removed then this only resets the current section's points
        // otherwise it creates new branch sections_ for each cluster and adds to the end of the sections_ list
        bifurcate(clusters);
        // update the max_radius, since the points have changed
        max_radius_ = radFromLength(sections_[sec_].max_distance_to_end);
      }
    }

    // we have split the ends, so we need to extract the set of nodes in a backwards manner
    if (extract_from_ends)
    {
      extractNodesFromEnds(nodes);
    }

    // estimate the section's tip (the centre of the cylinder of points)
    sections_[sec_].tip = calculateTipFromNodes(nodes);
    // get section's direction
    Eigen::Vector3d dir = par >= 0 ? (sections_[sec_].tip - sections_[par].tip).normalized() : Eigen::Vector3d(0, 0, 1);
    // shift to cylinder's centre
    sections_[sec_].tip += vectorToCylinderCentre(nodes, dir);
    // now find the segment radius
    sections_[sec_].radius = estimateCylinderRadius(nodes, dir);

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections_[sec_].ends.size() > 0)
    {
      addChildSection();
    }
  }  // end of loop. We now have created all of the BranchSections

  // Now calculate the section ids for all of the points, for the segmented cloud
  std::vector<int> section_ids(points_.size(), -1);
  calculateSectionIds(roots_list, section_ids, children);

  // debug draw to rviz the set of cylinders
  drawTrees(verbose);

  generateLocalSectionIds();

  Eigen::Vector3d min_bound(0, 0, 0), max_bound(0, 0, 0);
  // remove all sections with a root out of bounds, if we have gridded the cloud with an overlap
  if (params_->grid_width)
  {
    removeOutOfBoundSections(cloud, min_bound, max_bound);
  }

  std::vector<int> root_segs(cloud.ends.size(), -1);
  // now colour the ray cloud based on the segmentation
  segmentCloud(cloud, root_segs, section_ids);

  if (params_->grid_width)  // also remove rays from the segmented cloud
  {
    removeOutOfBoundRays(cloud, min_bound, max_bound, root_segs);
  }
}

/// calculate a branch radius from its length. We use allometric scaling based on a radius exponent and a
/// linear range for the final taper of the branch
/// See "Allometric patterns in Acer platanoides (Aceraceae) branches" for real world data
double Trees::radFromLength(double length)
{
  if (length / params_->linear_range)
  {
    return length / params_->length_to_radius;
  }
  const double radius_at_min = params_->linear_range / params_->length_to_radius;
  const double unscaled_rad = std::pow(radius_at_min, params_->radius_exponent);
  // length = scale_factor * rad^radius_exponent
  return std::pow(length * unscaled_rad / params_->linear_range, 1.0 / params_->radius_exponent);
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
    sections_[i].radius = radFromLength(sections_[i].max_distance_to_end);
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

Eigen::Vector3d Trees::getRootPosition()
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
// return the base point (average of roots)
void Trees::extractNodesAndEndsFromRoots(std::vector<int> &nodes, const Eigen::Vector3d &base,
                                         const std::vector<std::vector<int>> &children)
{
  const double thickness = params_->cylinder_length_to_width * max_radius_;
  const double thickness_sqr = thickness * thickness;
  const int par = sections_[sec_].parent;

  nodes = sections_[sec_].roots;
  // 2. find all the points (and the end points) for this section:
  for (unsigned int ijk = 0; ijk < nodes.size(); ijk++)
  {
    const int i = nodes[ijk];
    for (auto &child : children[i])
    {
      const double dist_sqr =
        par == -1 ? sqr(points_[child].pos[2] - base[2]) : (points_[child].pos - base).squaredNorm();
      if (dist_sqr < thickness_sqr)  // in same slot, so accumulate
      {
        nodes.push_back(child);  // so we recurse on this child too
      }
      else
      {
        sections_[sec_].ends.push_back(child);
      }
    }
  }
}

// find clusters of points from the root points up the shortest paths, up to the cylinder length
std::vector<std::vector<int>> Trees::findPointClusters(const Eigen::Vector3d &base, bool &points_removed)
{
  const double thickness = params_->cylinder_length_to_width * max_radius_;
  const int par = sections_[sec_].parent;

  std::vector<int> all_ends = sections_[sec_].ends;
  // 3. cluster end points to find if we have separate branches
  // first, get interpolated edge points_. i.e. interpolate between two connected points inside and outside
  // the branch section's cylinder, so that the set edge_pos are right on the top boundary of the section
  for (auto &j : all_ends)
  {
    const double dist1j = par == -1 ? points_[j].pos[2] - base[2] : (points_[j].pos - base).norm();
    const double dist0j =
      par == -1 ? points_[points_[j].parent].pos[2] - base[2] : (points_[points_[j].parent].pos - base).norm();
    const double blendj = (thickness - dist0j) / (dist1j - dist0j);
    points_[j].edge_pos = points_[points_[j].parent].pos * (1.0 - blendj) + points_[j].pos * blendj;
  }
  // convert to a structure that is better for the cluster function
  std::vector<Eigen::Vector3d> ps;
  std::vector<int> v_indices;
  for (auto &i : all_ends)
  {
    ps.push_back(points_[i].edge_pos);
    v_indices.push_back(i);
  }
  // cluster these end points based on two separation criteria (gap_ratio and span_ratio)
  std::vector<std::vector<int>> clusters;
  generateClusters(clusters, ps, params_->gap_ratio * max_radius_, params_->span_ratio * max_radius_);
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
void Trees::bifurcate(const std::vector<std::vector<int>> &clusters)
{
  const double thickness = params_->cylinder_length_to_width * max_radius_;
  const int par = sections_[sec_].parent;
  // find the maximum distance (to tip) for each cluster
  std::vector<double> max_distances(clusters.size());
  double maxmax = -1;
  int maxi = -1;
  for (size_t i = 0; i < clusters.size(); i++)
  {
    max_distances[i] = 0;
    for (auto &end : clusters[i])
    {
      max_distances[i] = std::max(max_distances[i], points_[end].distance_to_end);
    }
    if (max_distances[i] > maxmax)
    {
      maxmax = max_distances[i];
      maxi = static_cast<int>(i);
    }
  }

  // set the current branch section to the cluster with the
  // maximum distance to tip (i.e. the longest branch).
  if (maxi == -1)
  {
    std::cout << "error: bad maxi" << std::endl;
  }
  sections_[sec_].ends = clusters[maxi];
  sections_[sec_].max_distance_to_end = max_distances[maxi] + thickness;

  // for all other clusters, add new sections to the list...
  for (size_t i = 0; i < clusters.size(); i++)
  {
    if (static_cast<int>(i) == maxi)
    {
      continue;
    }
    BranchSection new_node = sections_[sec_];
    new_node.max_distance_to_end = max_distances[i] + thickness;
    const double maxrad = radFromLength(new_node.max_distance_to_end);
    if (maxrad > 0.5 * params_->min_diameter)  // but only add if they are large enough
    {
      // we only specify the end points at this stage. They will therefore enter the
      // extract_from_ends block below when the sec_ for loop gets to their section
      new_node.ends = clusters[i];
      if (par != -1)
      {
        sections_[par].children.push_back(static_cast<int>(sections_.size()));
      }
      sections_.push_back(new_node);
    }
  }
}

// find the nodes between the section end points and the section root points
void Trees::extractNodesFromEnds(std::vector<int> &nodes)
{
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

// use the nodes to estimate a tip location (which is the mean of the node points)
Eigen::Vector3d Trees::calculateTipFromNodes(const std::vector<int> &nodes)
{
  Eigen::Vector3d tip(0, 0, 0);
  if (nodes.empty())
  {
    std::cout << "error: there shouldn't be empty nodes at this point in the processing" << std::endl;
  }
  // get the points in this segment
  auto &list = nodes.empty() ? (sections_[sec_].ends.empty() ? sections_[sec_].roots : sections_[sec_].ends) : nodes;

  for (auto &i : list)
  {
    tip += points_[i].pos;
  }
  if (list.size() > 0)
  {
    tip /= static_cast<double>(list.size());
  }
  return tip;
}

Eigen::Vector3d Trees::vectorToCylinderCentre(const std::vector<int> &nodes, const Eigen::Vector3d &dir)
{
#define REAL_CENTROID  // finds a new centroid that is robust to branches scanned from a single side
#if defined REAL_CENTROID
  const int par = sections_[sec_].parent;
  Eigen::Vector3d mean_p(0, 0, 0);
  std::vector<Eigen::Vector3d> ps;
  const Eigen::Vector3d vec(1, 2, 3);
  // obtain two orthogonal planes to the section's direction vector
  const Eigen::Vector3d ax1 = dir.cross(vec).normalized();
  const Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
  double n = 0;
  // iterate over the nodes and the ends points_...
  const int start_ii = par >= 0 || sections_[sec_].ends.empty() ? 0 : 1;
  for (int ii = start_ii; ii < 2; ii++)
  {
    const std::vector<int> &node_list = ii == 0 ? nodes : sections_[sec_].ends;
    // project points into a local space parabola
    for (auto &i : node_list)
    {
      const Eigen::Vector3d pos = points_[i].pos - sections_[sec_].tip;
      const Eigen::Vector2d offset(ax1.dot(pos), ax2.dot(pos));
      const Eigen::Vector2d xy = offset / sections_[sec_].radius;
      const double l2 = xy.squaredNorm();
      const Eigen::Vector3d point(xy[0], xy[1], 0.5 * l2);  // a paraboloid that has gradient 1 at 1
      ps.push_back(point);
      mean_p += point;
      n++;
    }
  }
  mean_p /= n;
  if (n > 5)  // assuming there are sufficient points for a resonable guess
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
      const double A = (plane.xz * plane.y2 - plane.yz * plane.xy) / (plane.x2 * plane.y2 - plane.xy * plane.xy);
      const double B = (plane.yz - A * plane.xy) / plane.y2;

      Eigen::Vector2d shift(A, B);
      const double shift2 = shift.squaredNorm();
      if (par >= 0 && shift2 > 1.0)  // don't shift more than one radius each iteration, for safety
      {
        shift /= std::sqrt(shift2);
      }

      // apply the plane parameters as a world-space shift in the branch section tip position
      return (ax1 * shift[0] + ax2 * shift[1]) * sections_[sec_].radius;
    }
  }
  else if (n > 0)  // if only a few points then use the mean instead
  {
    return (ax1 * mean_p[0] + ax2 * mean_p[1]) * sections_[sec_].radius;
  }
#endif
  return Eigen::Vector3d(0, 0, 0);
}

double Trees::estimateCylinderRadius(const std::vector<int> &nodes, const Eigen::Vector3d &dir)
{
  int par = sections_[sec_].parent;
  if (par == -1)  // if this is the root segment
  {
    double n = 0, rad = 0;
    // then get the mean radius
    const auto &list = sections_[sec_].ends.empty() ? nodes : sections_[sec_].ends;
    for (auto &id : list)
    {
      const Eigen::Vector3d offset = points_[id].pos - sections_[sec_].tip;
      rad += (offset - dir * offset.dot(dir)).norm();
      n++;
    }
    const double radius = list.size() < 2 ? max_radius_ : rad / n;
    return std::min(radius, max_radius_);
  }
  // for non-root segments

  // use the parent radius as a prior with a weight of 4 points
  // this avoids spurious radius estimations when the number of points is
  // small
  double n = 4.0;
  double rad = n * sections_[par].radius;
  for (auto &node : nodes)
  {
    const Eigen::Vector3d offset = points_[node].pos - sections_[sec_].tip;
    rad += (offset - dir * offset.dot(dir)).norm();
    n++;
  }
  return std::min(std::min(rad / n, max_radius_), sections_[par].radius);
}

// add a child section to continue reconstructing the tree segments
void Trees::addChildSection()
{
  BranchSection new_node;
  new_node.parent = static_cast<int>(sec_);
  new_node.roots = sections_[sec_].ends;
  new_node.max_distance_to_end = 0.0;
  for (auto &root : new_node.roots)
  {
    new_node.max_distance_to_end = std::max(new_node.max_distance_to_end, points_[root].distance_to_end);
  }
  max_radius_ = radFromLength(new_node.max_distance_to_end);
  // constrain each new node's radius to be no larger than its parent radius. This is a reasonable constraint.
  new_node.radius = std::min(sections_[sec_].radius, max_radius_);
  if (new_node.radius > 0.5 * params_->min_diameter)  // if it is the first node, then we need a second node
  {
    sections_[sec_].children.push_back(static_cast<int>(sections_.size()));

    new_node.tip.setZero();
    for (auto &end : new_node.roots)
    {
      new_node.tip += points_[end].pos;
    }
    if (new_node.roots.size() > 0)
    {
      new_node.tip /= static_cast<double>(new_node.roots.size());
    }
    sections_.push_back(new_node);
  }
}

// calcualte what section every point belongs to
void Trees::calculateSectionIds(const std::vector<std::vector<int>> &roots_list, std::vector<int> &section_ids,
                                const std::vector<std::vector<int>> &children)
{
  // first we initialise to one branch section per root
  for (size_t i = 0; i < roots_list.size(); i++)
  {
    for (auto &root : roots_list[i])
    {
      section_ids[root] = static_cast<int>(i);
    }
  }
  for (sec_ = 0; sec_ < sections_.size(); sec_++)
  {
    std::vector<int> nodes;
    if (sections_[sec_].ends.size() > 0)
    {
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

void Trees::drawTrees(bool verbose)
{
  if (verbose)
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<double> radii;
    for (auto &tree_node : sections_)
    {
      if (tree_node.tip == Eigen::Vector3d(0, 0, 0))
      {
        continue;
      }
      if (tree_node.parent >= 0)
      {
        if ((sections_[tree_node.parent].tip - tree_node.tip).norm() < 0.001)
        {
          continue;
        }
        starts.push_back(sections_[tree_node.parent].tip);
        ends.push_back(tree_node.tip);
        radii.push_back(tree_node.radius);
      }
    }
    DebugDraw::instance()->drawCylinders(starts, ends, radii, 0);
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

void Trees::removeOutOfBoundSections(const Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound)
{
  const double width = params_->grid_width;
  cloud.calcBounds(&min_bound, &max_bound);
  const Eigen::Vector3d mid = (min_bound + max_bound) / 2.0;
  const Eigen::Vector2i inds(std::round(mid[0] / width), std::round(mid[1] / width));
  min_bound[0] = width * (static_cast<double>(inds[0]) - 0.5);
  min_bound[1] = width * (static_cast<double>(inds[1]) - 0.5);
  max_bound[0] = width * (static_cast<double>(inds[0]) + 0.5);
  max_bound[1] = width * (static_cast<double>(inds[1]) + 0.5);
  std::cout << "min bound: " << min_bound.transpose() << ", max bound: " << max_bound.transpose() << std::endl;

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

      if (!params_->segment_branches)
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
}

// remove rays from the ray cloud where the end points are out of bounds
void Trees::removeOutOfBoundRays(Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound,
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
bool Trees::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << "# tree structure file:" << std::endl;
  ofs << "x,y,z,radius,parent_id,section_id" << std::endl;  // simple format
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent >= 0 || section.children.empty())  // not a root section, so move on
    {
      continue;
    }
    ofs << section.tip[0] << "," << section.tip[1] << "," << section.tip[2] << "," << section.radius << ",-1," << sec;

    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      const BranchSection &node = sections_[children[c]];
      ofs << ", " << node.tip[0] << "," << node.tip[1] << "," << node.tip[2] << "," << node.radius << ","
          << sections_[node.parent].id;
      ofs << ", " << children[c];
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
