// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raysegment.h"
#include <nabo/nabo.h>
#include "rayterrain.h"
#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>

namespace ray
{
/// nodes of priority queue used in shortest path algorithm
struct QueueNode
{
//  QueueNode() {}
  QueueNode(double distance_to_ground, double score, double radius, int root, int index)
    : distance_to_ground(distance_to_ground)
    , score(score)
    , radius(radius)
    , root(root)
    , id(index)
  {}

  double distance_to_ground; // path distance to the ground
  double score;              // score is the modified edge length metric being minimised
  double radius;             // radius of the tree base, this acts as a score scale coefficient
  int root;                  // index of the root of the path
  int id;                    // index into the points_ array for this node
};

class QueueNodeComparator
{
public:
  bool operator()(const QueueNode &p1, const QueueNode &p2) { return p1.score > p2.score; }
};

/// Connect the supplied set of points @c points according to the shortest path to the ground, by filling in their
/// parent indices
/// @c distance_limit maximum distance between points that can be connected
/// @c gravity_factor controls how far laterally the shortest paths can travel
/// @c closest_node a priority queue
/// @c path_labels optional per-point tree labels (-1 for unlabelled). When supplied, paths are not allowed to
/// connect points of two different labels, and unlabelled points adopt the label of the path that reaches them
void connectPointsShortestPath(
  std::vector<Vertex> &points,
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node, double distance_limit,
  double gravity_factor, std::vector<int> *path_labels)
{
  // 1. get nearest neighbours
  const int search_size = std::min(20, static_cast<int>(points.size()) - 1);
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++)
  {
    points_p.col(i) = points[i].pos;
  }
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, distance_limit);

  // 2. climb up from lowest points, this part is based on Djikstra's algorithm
  while (!closest_node.empty())
  {
    QueueNode node = closest_node.top();
    closest_node.pop();
    if (!points[node.id].visited)
    {
      // for each unvisited point, look at its nearest neighbours
      for (int i = 0; i < search_size && indices(i, node.id) != Nabo::NNSearchD::InvalidIndex; i++)
      {
        const int child = indices(i, node.id);
        if (path_labels)
        {
          // paths are not allowed to cross from one labelled tree into another
          const int node_label = (*path_labels)[node.id];
          const int child_label = (*path_labels)[child];
          if (node_label != -1 && child_label != -1 && node_label != child_label)
          {
            continue;
          }
        }
        const double dist2 = dists2(i, node.id);  // square distance to neighbour
        const double dist = std::sqrt(dist2);
        double new_score = 0;
        const Eigen::Vector3d dif = (points[child].pos - points[node.id].pos).normalized();
        Eigen::Vector3d dir(0, 0, 1);
        // estimate direction of path from parent of parent if possible
        const int ppar = points[node.id].parent;
        if (ppar != -1)
        {
          if (points[ppar].parent != -1)  // this is a bit smoother than...
          {
            dir = (points[node.id].pos - points[points[ppar].parent].pos).normalized();
          }
          else  // ..just this
          {
            dir = (points[node.id].pos - points[ppar].pos).normalized();
          }
        }
        const double d = std::max(0.001, dif.dot(dir));
        // we are looking for a minimum score, so large distances are bad, but new points in line with the
        // path direction are good
        double score = dist2 / (d * d);

        if (gravity_factor > 0.0)  // penalise paths that are hard to hold up against gravity (lateral direction)
        {
          Eigen::Vector3d to_node = points[node.id].pos - points[node.root].pos;
          to_node[2] = 0.0;
          const double lateral_sqr = to_node.squaredNorm();
          const double gravity_scale =
            1.0 + gravity_factor * lateral_sqr;  // the squaring means gravity plays little role for normal trees,
                                                 // kicking in stronger on outlier lateral ones
          score *= gravity_scale;
        }

        // scale score according to size of each tree, this prevents small trees from
        // capturing the branches of larger trees
        score /= node.radius;

        new_score = node.score + score;
        if (new_score < points[child].score)
        {
          points[child].score = new_score;
          // we also maintain the distance to ground value
          points[child].distance_to_ground = node.distance_to_ground + dist;
          points[child].parent = node.id;
          points[child].root = node.root;
          if (path_labels && (*path_labels)[child] == -1)
          {
            // an unlabelled point belongs to the tree that its path comes from
            (*path_labels)[child] = (*path_labels)[node.id];
          }
          closest_node.push(
            QueueNode(points[child].distance_to_ground, points[child].score, node.radius, node.root, child));
        }
      }
      points[node.id].visited = true;
    }
  }
}

/// Converts a ray cloud to a set of points @c points connected by the shortest path to the ground @c mesh
/// the returned vector of index sets provides the root points for each separated tree
std::vector<std::vector<int>> getRootsAndSegment(std::vector<Vertex> &points, const Cloud &cloud, const Mesh &mesh,
                                                 double max_diameter, double distance_limit, double height_min,
                                                 double gravity_factor, const std::vector<int> *point_labels)
{
  if (point_labels && point_labels->size() != cloud.ends.size())
  {
    std::cerr << "Error: the number of point labels does not match the number of rays in the cloud" << std::endl;
    return std::vector<std::vector<int>>();
  }
  // first fill in the basic attributes of the points structure
  points.reserve(cloud.ends.size());
  std::vector<int> labels;  // the supplied label per point, or -1 where unlabelled
  int num_labels = 0;
  if (point_labels)
  {
    labels.reserve(cloud.ends.size());
  }
  // maps the supplied label values onto a contiguous set of tree indices, in ascending label order so that
  // the output trees keep the order (and so the IDs) of the input labels
  std::map<int, int> label_ids;
  if (point_labels)
  {
    for (unsigned int i = 0; i < cloud.ends.size(); i++)
    {
      if (cloud.rayBounded(i) && (*point_labels)[i] >= 0)
      {
        label_ids[(*point_labels)[i]] = 0;
      }
    }
    for (auto &label_id : label_ids)
    {
      label_id.second = num_labels++;
    }
  }
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    if (cloud.rayBounded(i))
    {
      points.push_back(Vertex(cloud.ends[i], cloud.starts[i]));
      if (point_labels)
      {
        // negative labels are unlabelled, so these points attach to whichever tree's path reaches them
        const int label = (*point_labels)[i];
        labels.push_back(label < 0 ? -1 : label_ids[label]);
      }
    }
  }
  const int num_cloud_points = static_cast<int>(points.size());
  if (point_labels)
  {
    std::cout << "using " << num_labels << " pre-labelled trees from the input cloud" << std::endl;
  }

  const double pixel_width = max_diameter;
  Eigen::Vector3d box_min, box_max;
  cloud.calcBounds(&box_min, &box_max);
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;

  // also add points for every vertex on the ground mesh.
  const int roots_start = static_cast<int>(points.size());
  for (auto &vert : mesh.vertices())
  {
    if (vert[0] >= box_min[0] && vert[1] >= box_min[1] &&
        vert[0] <= box_max[0] && vert[1] <= box_max[1])
    {
      points.push_back(Vertex(vert - Eigen::Vector3d(0,0,0.01), vert + Eigen::Vector3d(0,0,0.01))); // make small vertical ray
    }
  }
  // the ground mesh vertices are unlabelled, they are the seeds that the labelled trees grow from
  std::vector<int> path_labels;
  if (point_labels)
  {
    path_labels = labels;
    path_labels.resize(points.size(), -1);
  }
  // convert the ground mesh to an easy look-up height field
  Eigen::ArrayXXd lowfield;
  mesh.toHeightField(lowfield, box_min, box_max, pixel_width);

  // set heightfield as the height of the canopy above the ground
  Eigen::ArrayXXd heightfield =
    Eigen::ArrayXXd::Constant(static_cast<int>(lowfield.rows()), static_cast<int>(lowfield.cols()), std::numeric_limits<double>::lowest());
  for (const auto &point : points)
  {
    Eigen::Vector3i index = ((point.pos - box_min) / pixel_width).cast<int>();
    heightfield(index[0], index[1]) = std::max(heightfield(index[0], index[1]), point.pos[2]);
  }
  // make heightfield relative to the ground
  for (int i = 0; i < heightfield.rows(); i++)
  {
    for (int j = 0; j < heightfield.cols(); j++)
    {
      heightfield(i, j) = std::max(1e-10, heightfield(i, j) - lowfield(i, j));
    }
  }

  // create an initial priority queue node for each root point (mesh vertex) using the
  // observed height as a scaling parameter
  for (int ind = roots_start; ind < static_cast<int>(points.size()); ind++)
  {
    points[ind].distance_to_ground = 0.0;
    points[ind].score = 0.0;
    points[ind].root = ind;
    const Eigen::Vector3i index = ((points[ind].pos - box_min) / pixel_width).cast<int>();
    closest_node.push(QueueNode(0, 0, heightfield(index[0], index[1]), ind, ind));
  }

  // perform Djikstra's shortest path to ground algorithm to fill in the parent indices in 'points'
  connectPointsShortestPath(points, closest_node, distance_limit, gravity_factor,
                            point_labels ? &path_labels : nullptr);

  // when labels have been supplied we don't need to segment the trees, we just group the root points
  // according to the label of the points that grow from them, generating one tree per label
  if (point_labels)
  {
    // each labelled point votes for its root point, so that a root point shared by several
    // labels is given to the label that most of its points belong to
    std::vector<std::unordered_map<int, int>> votes(points.size() - roots_start);
    std::vector<double> max_heights(num_labels, 0.0);
    for (int i = 0; i < num_cloud_points; i++)
    {
      const int label = labels[i];
      if (label == -1)
      {
        continue;
      }
      const Eigen::Vector3i index = ((points[i].pos - box_min) / pixel_width).cast<int>();
      max_heights[label] = std::max(max_heights[label], points[i].pos[2] - lowfield(index[0], index[1]));
      if (points[i].root >= roots_start)
      {
        votes[points[i].root - roots_start][label]++;
      }
    }
    std::vector<std::vector<int>> roots_set(num_labels);
    std::vector<int> best_root(num_labels, -1);       // the most voted for root point per label
    std::vector<int> best_root_votes(num_labels, 0);
    for (int i = roots_start; i < static_cast<int>(points.size()); i++)
    {
      int best_label = -1;
      int most_votes = 0;
      for (const auto &vote : votes[i - roots_start])
      {
        if (vote.second > most_votes)
        {
          most_votes = vote.second;
          best_label = vote.first;
        }
      }
      if (best_label != -1)
      {
        roots_set[best_label].push_back(i);
        if (most_votes > best_root_votes[best_label])
        {
          best_root_votes[best_label] = most_votes;
          best_root[best_label] = i;
        }
      }
    }
    // only keep the root points close to the most voted for one, so that a single path descending to a
    // distant patch of ground cannot move the base of the tree
    const double max_root_radius = 2.0 * pixel_width;
    for (int label = 0; label < num_labels; label++)
    {
      if (best_root[label] == -1)
      {
        continue;
      }
      std::vector<int> close_roots;
      for (auto &root : roots_set[label])
      {
        Eigen::Vector3d dif = points[root].pos - points[best_root[label]].pos;
        dif[2] = 0.0;
        if (dif.norm() <= max_root_radius)
        {
          close_roots.push_back(root);
        }
      }
      roots_set[label] = close_roots;
    }
    // a label whose points all descend through another label's root points has no root of its own,
    // so we use its entry points (the points where its paths cross into the label) instead
    std::vector<bool> needs_entry_points(num_labels);
    for (int label = 0; label < num_labels; label++)
    {
      needs_entry_points[label] = roots_set[label].empty();
    }
    for (int i = 0; i < num_cloud_points; i++)
    {
      const int label = labels[i];
      if (label == -1 || !needs_entry_points[label])
      {
        continue;
      }
      const int parent = points[i].parent;
      if (parent == -1 || path_labels[parent] != label)
      {
        roots_set[label].push_back(i);
      }
    }
    // lastly remove the labels that are too small to count as a tree
    std::vector<std::vector<int>> result;
    int num_short = 0;
    for (int label = 0; label < num_labels; label++)
    {
      if (roots_set[label].empty())
      {
        continue;
      }
      if (max_heights[label] < height_min)
      {
        num_short++;
        continue;
      }
      result.push_back(roots_set[label]);
    }
    if (num_short > 0)
    {
      std::cout << "removed " << num_short << " labelled trees shorter than the " << height_min
                << " m minimum height" << std::endl;
    }
    return result;
  }

  // next we want to segment the paths into separate trees. To do this we find the number of points and
  // the maximum height of points that come from each cell index
  Eigen::ArrayXXi counts =
    Eigen::ArrayXXi::Constant(static_cast<int>(heightfield.rows()), static_cast<int>(heightfield.cols()), 0);
  Eigen::ArrayXXd heights =
    Eigen::ArrayXXd::Constant(static_cast<int>(heightfield.rows()), static_cast<int>(heightfield.cols()), 0);
  for (const auto &point : points)
  {
    if (point.root == -1)
    {
      continue;
    }
    const Eigen::Vector3i index = ((points[point.root].pos - box_min) / pixel_width).cast<int>();
    counts(index[0], index[1])++;
    heights(index[0], index[1]) = std::max(heights(index[0], index[1]), point.pos[2] - lowfield(index[0], index[1]));
  }

  // in order to avoid boundary artefacts, we create a 2x2 summed array:
  Eigen::ArrayXXi sums = Eigen::ArrayXXi::Constant(static_cast<int>(counts.rows()), static_cast<int>(counts.cols()), 0);
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      const int i2 = std::min(i + 1, static_cast<int>(sums.rows()) - 1);
      const int j2 = std::min(j + 1, static_cast<int>(sums.cols()) - 1);
      sums(i, j) = counts(i, j) + counts(i, j2) + counts(i2, j) + counts(i2, j2);
    }
  }
  // now find the best 2x2 sum for each cell:
  std::vector<Eigen::Vector2i> bests(static_cast<int>(counts.rows()) * static_cast<int>(counts.cols()));
  for (int x = 0; x < static_cast<int>(sums.rows()); x++)
  {
    for (int y = 0; y < static_cast<int>(sums.cols()); y++)
    {
      Eigen::Vector2i best_index(-1, -1);
      int largest_sum = -1;
      for (int i = std::max(0, x - 1); i <= x; i++)
      {
        for (int j = std::max(0, y - 1); j <= y; j++)
        {
          if (sums(i, j) > largest_sum)
          {
            largest_sum = sums(i, j);
            best_index = Eigen::Vector2i(i, j);
          }
        }
      }
      bests[x + sums.rows() * y] = best_index;
    }
  }
  // next we need to find the highest point for each cell....
  Eigen::ArrayXXd max_heights =
    Eigen::ArrayXXd::Constant(static_cast<int>(counts.rows()), static_cast<int>(counts.cols()), 0);
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      Eigen::Vector2i best_index = bests[i + sums.rows() * j];
      double max_height = 0.0;
      for (int x = best_index[0]; x < std::min(best_index[0] + 2, static_cast<int>(sums.rows())); x++)
      {
        for (int y = best_index[1]; y < std::min(best_index[1] + 2, static_cast<int>(sums.cols())); y++)
        {
          if (bests[x + static_cast<int>(sums.rows()) * y] == best_index)
          {
            max_height = std::max(max_height, heights(x, y));
          }
        }
      }
      max_heights(best_index[0], best_index[1]) = max_height;
    }
  }

  // now that we have a max height for each cell, we can fill in a list of the root
  // points for each cell
  std::vector<std::vector<int>> roots_lists(sums.rows() * sums.cols());
  for (int i = roots_start; i < static_cast<int>(points.size()); i++)
  {
    const Eigen::Vector3i index = ((points[i].pos - box_min) / pixel_width).cast<int>();
    const Eigen::Vector2i best_index = bests[index[0] + static_cast<int>(sums.rows()) * index[1]];
    const double max_height = max_heights(best_index[0], best_index[1]);
    if (max_height >= height_min)
    {
      const int id = best_index[0] + static_cast<int>(sums.rows()) * best_index[1];
      roots_lists[id].push_back(i);
    }
  }

  // convert this into a contiguous form, which represents the set of
  // root points for each tree
  std::vector<std::vector<int>> roots_set;
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      auto &roots = roots_lists[i + static_cast<int>(sums.rows()) * j];
      if (roots.size() > 0)
      {
        roots_set.push_back(roots);
      }
    }
  }

  return roots_set;
}

}  // namespace ray
