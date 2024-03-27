// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrunks.h"
#include <nabo/nabo.h>
#include <queue>
#include "../raycuboid.h"
#include "../raygrid.h"
#include "../rayply.h"
#include "raygrid2d.h"

namespace ray
{
namespace
{
/// debug draw the trunks, either as a set of lines to rviz, or as a point cloud.
void drawTrunks(const std::vector<Trunk> &trunks, const Eigen::Vector3d &offset, const std::vector<Eigen::Vector3d> *closest_approach_points = nullptr,
                const std::vector<Eigen::Vector3d> *pass_through_points = nullptr)
{
  std::vector<Eigen::Vector3d> cloud_points;
  std::vector<double> times;
  std::vector<RGBA> colours;
  RGBA colour;
  colour.red = colour.green = colour.blue = 0;
  colour.alpha = 255;

  // add the optional additional points to the ray cloud
  if (closest_approach_points)
  {
    for (auto &point : *closest_approach_points)
    {
      cloud_points.push_back(point + offset);
      times.push_back(0.0);
      colours.push_back(colour);
    }
  }
  colour.red = colour.green = 255;
  if (pass_through_points)
  {
    for (auto &point : *pass_through_points)
    {
      cloud_points.push_back(point + offset);
      times.push_back(0.0);
      colours.push_back(colour);
    }
  }
  // Now draw each of the trunks as a cylinder composed of a set of points
  for (auto const &trunk : trunks)
  {
    if (!trunk.active)
    {
      continue;
    }
    // grey-scale colour where black to white represents a score of 0 to twice the minimum score.
    // this range allows results above and below the minimum score to be clearly discerned. 
    colour.red = colour.green = colour.blue = (uint8_t)std::min(255.0 * trunk.score / (2.0 * minimum_score), 255.0);
    // in order to distinguish results below the minimum score, we keep only the red component. 
    // red, as is typical, represents failure (to meet the criteria of a cylindrical trunk)
    if (trunk.score < minimum_score)
    {
      colour.green = colour.blue = 0;
    }

    const Eigen::Vector3d side1 = trunk.dir.cross(Eigen::Vector3d(1, 2, 3)).normalized();
    const Eigen::Vector3d side2 = side1.cross(trunk.dir);
    const double rad = trunk.radius;
    // the point-cloud cylinder is composed of...
    for (double z = -0.5; z < 0.5; z += 0.15)  // multiple heights
    {
      for (double ang = 0.0; ang < 2.0 * kPi; ang += 0.6)  // and multiple angles
      {
        const Eigen::Vector3d pos =
          trunk.centre + trunk.dir * trunk.length * z + side1 * std::sin(ang) * rad + side2 * std::cos(ang) * rad;
        cloud_points.push_back(pos + offset);
        times.push_back(0.0);
        colours.push_back(colour);
      }
    }
  }
  writePlyPointCloud("trunks_verbose.ply", cloud_points, times, colours);
}

/// Used in a priority queue for identifying close neighbouring trunks
struct QueueNode
{
  QueueNode(double score, int index)
    : score(score)
    , id(index)
  {}

  double score;
  int id;
};

class QueueNodeComparator
{
public:
  bool operator()(const QueueNode &p1, const QueueNode &p2) { return p1.score > p2.score; }
};
}  // namespace

// A map of voxels to integers, such as for a count value per voxel
class IntegerVoxels
{
public:
  /// Construct from a given @c width and vector @c offset
  IntegerVoxels(double width, const Eigen::Vector3d &offset)
    : voxel_width(width)
    , offset(offset)
  {}

  /// Obtain the 3D index that overlaps any 3D position @c pos
  inline Eigen::Vector3i getIndex(const Eigen::Vector3d &pos)
  {
    const Eigen::Vector3d ind = (pos - offset) / voxel_width;
    return Eigen::Vector3i(int(std::floor(ind[0])), int(std::floor(ind[1])), int(std::floor(ind[2])));
  }

  /// For a given 3D position @c pos, increment the value of the overlapping voxel
  inline void increment(const Eigen::Vector3d &pos) { increment(getIndex(pos)); }

  /// increment the value of the voxel at the given 3D @c index
  inline void increment(const Eigen::Vector3i &index)
  {
    auto it = voxel_map.find(index);
    if (it == voxel_map.end())
    {
      voxel_map[index] = 1;
    }
    else
    {
      it->second++;
    }
  }

  /// get the integer value of the voxel at the specified 3D @c index
  inline int get(const Eigen::Vector3i &index) const
  {
    auto it = voxel_map.find(index);
    if (it == voxel_map.end())
    {
      return 0;
    }
    else
    {
      return it->second;
    }
  }

  /// apply the specified function @c func to every voxel in the grid
  void forEach(std::function<void(const struct IntegerVoxels &voxels, double width, const Eigen::Vector3d &offset,
                                  const Eigen::Vector3i &index, int count)>
                 func)
  {
    for (auto &voxel : voxel_map)
    {
      func(*this, voxel_width, offset, voxel.first, voxel.second);
    }
  }
private:
  std::map<Eigen::Vector3i, int, Vector3iLess> voxel_map;
  double voxel_width;
  Eigen::Vector3d offset;
};

/// This method generates a set of initial trunk candidates at four different scales, by gridding the 
/// space and initialising a vertically placed trunk centred in each voxel that contains sufficient points
/// in its neighbourhoood
void initialiseTrunks(std::vector<Trunk> &trunks, const Cloud &cloud, const Eigen::Vector3d &min_bound,
                      double voxel_width)
{
  // we generate a set of candidate trunks at multiple resolutions
  const Eigen::Vector3d half_voxel(0.5 * voxel_width, 0.5 * voxel_width, 0.5 * voxel_width);
  // normal size, with offset grid
  IntegerVoxels voxels1(voxel_width, min_bound), voxels2(voxel_width, min_bound + half_voxel);
  // double size with offset grid
  IntegerVoxels voxels3(2.0 * voxel_width, min_bound), voxels4(2.0 * voxel_width, min_bound + 2.0 * half_voxel);
  // this is used to count the number of points in each voxel
  IntegerVoxels *voxels[4] = { &voxels1, &voxels2, &voxels3, &voxels4 };
  Eigen::Vector3i ind2;
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    for (int j = 0; j < 4; j++)
    {
      voxels[j]->increment(cloud.ends[i]);
    }
  }
  // add a new trunk for each voxel
  auto addTrunk = [&trunks](const struct IntegerVoxels &voxels, double voxel_width, const Eigen::Vector3d &offset,
                            const Eigen::Vector3i &index, int count) {
    const Eigen::Vector3d centre = (index.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width + offset;
    if (count < 4)  // not enough points to add a candidate trunk
    {
      return;
    }
    // now check that there are points above and below the voxel in question. To make sure that there could be
    // a tall enough trunk here in the data
    const int minpoints = 4;
    const bool above1 = voxels.get(index + Eigen::Vector3i(0, 0, 1)) >= minpoints;
    const bool above2 = voxels.get(index + Eigen::Vector3i(0, 0, 2)) >= minpoints;
    const bool below1 = voxels.get(index - Eigen::Vector3i(0, 0, 1)) >= minpoints;
    const bool below2 = voxels.get(index - Eigen::Vector3i(0, 0, 2)) >= minpoints;
    const bool tall_enough = (above1 && below1) || (above1 && above2) || (below1 && below2);
    if (!tall_enough)
    {
      return;
    }
    // ok this is good enough for a candidate, so add it in
    Trunk trunk;
    const double diameter = voxel_width / std::sqrt(2.0);
    trunk.centre = centre;
    trunk.radius = diameter / 2.0;
    trunk.length = diameter * trunk_height_to_width;
    trunk.score = 0;
    trunk.dir = Eigen::Vector3d(0, 0, 1);
    trunks.push_back(trunk);
  };
  for (int j = 0; j < 4; j++)
  {
    // now run the addTrunk function for each voxel
    voxels[j]->forEach(addTrunk);
  }
}

// get rid of trunks that overlap existing trunks
void removeOverlappingTrunks(std::vector<Trunk> &best_trunks_)
{
  // Brute force approach at the moment!
  for (size_t i = 0; i < best_trunks_.size(); i++)
  {
    if (!(i % 1000))
    {
      std::cout << i << " / " << best_trunks_.size() << std::endl;
    }
    Trunk &trunk = best_trunks_[i];
    if (!trunk.active)
    {
      continue;
    }

    // ideally we would get the volume of overlap of two cylinders, however
    // since this is complicated, we approximate by taking a multiple points in the
    // first cylinder and counting how many overlap the second cylinder
    const Eigen::Vector3d ax1 = Eigen::Vector3d(1, 2, 3).cross(trunk.dir).normalized();
    const Eigen::Vector3d ax2 = trunk.dir.cross(ax1);
    const int num = 5;
    const double s = 0.8;
    // 5 points in a cross
    const double xs[num] = { 0, s, 0, -s, 0 };
    const double ys[num] = { 0, 0, s, 0, -s };
    // x 5 heights along the cylinder = 25 points
    const double zs[num] = { -0.5 * s, -0.25 * s, 0, 0.25 * s, 0.5 * s };

    // for coarse intersection
    const Eigen::Vector3d base = trunk.centre - 0.5 * trunk.length * trunk.dir;
    const Eigen::Vector3d top = trunk.centre + 0.5 * trunk.length * trunk.dir;
    const Eigen::Vector3d rad(trunk.radius, trunk.radius, trunk.radius);
    const Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

    for (size_t j = 0; j < best_trunks_.size(); j++)  // brute force
    {
      if (i == j)
      {
        continue;
      }
      Trunk &cylinder = best_trunks_[j];
      if (!cylinder.active)
      {
        continue;
      }
      // first, bounding box intersection test
      const Eigen::Vector3d base2 = cylinder.centre - 0.5 * cylinder.length * cylinder.dir;
      const Eigen::Vector3d top2 = cylinder.centre + 0.5 * cylinder.length * cylinder.dir;
      const Eigen::Vector3d rad2(cylinder.radius, cylinder.radius, cylinder.radius);
      const Cuboid cuboid2(minVector(base2, top2) - rad2, maxVector(base2, top2) + rad2);
      if (!cuboid.overlaps(cuboid2))  // broadphase exclusion
      {
        continue;
      }

      // now count the number of intersections
      int num_inside = 0;
      std::vector<Eigen::Vector3d> ps;
      for (int k = 0; k < num; k++)
      {
        for (int l = 0; l < num; l++)
        {
          Eigen::Vector3d pos =
            trunk.centre + trunk.dir * zs[k] * trunk.length + (ax1 * xs[l] + ax2 * ys[l]) * trunk.radius;
          pos -= cylinder.centre;
          // is pos inside cylinder j?
          const double d = pos.dot(cylinder.dir);
          if (d > cylinder.length * 0.5 || d < -cylinder.length * 0.5)
          {
            continue;
          }
          pos -= cylinder.dir * d;
          if (pos.squaredNorm() < sqr(cylinder.radius))
          {
            num_inside++;
          }
        }
      }

      const double inside_ratio = static_cast<double>(num_inside) / static_cast<double>(num * num);
      if (inside_ratio > 0.4)
      {
        const double vol_trunk = sqr(trunk.radius) * trunk.length;
        const double vol_cylinder = sqr(cylinder.radius) * cylinder.length;
        // remove the smaller one
        if (vol_trunk < vol_cylinder)
        {
          trunk.active = false;
        }
        else
        {
          cylinder.active = false;
        }
        break;
      }
    }
  }
}

// trunk identification in ray cloud
Trunks::Trunks(const Cloud &cloud, const Eigen::Vector3d &offset, double midRadius, bool verbose, bool remove_permeable_trunks)
{
  // The method is iterative, starting with a large set of trunk candidates, it
  // iteratively adjusts their pose and size to better approximate the neighbourhood of points,
  // filtering out poor candidates on the fly.
  // after the iteration loop has finished, the trunks are further filtered to remove overlapping
  // candidates, and the result is saved to best_trunks_

  // spacing helps us choose a grid resolution that is appropriate to the cloud
  const double spacing = cloud.estimatePointSpacing();
  if (verbose)
  {
    std::cout << "av radius: " << midRadius << ", estimated point spacing: " << spacing
              << ", minimum score: " << minimum_score << std::endl;
  }
  const Eigen::Vector3d min_bound = cloud.calcMinBound();
  const Eigen::Vector3d max_bound = cloud.calcMaxBound();
  std::cout << "cloud from: " << min_bound.transpose() << " to: " << max_bound.transpose() << std::endl;

  std::vector<Trunk> trunks;

  // 1. voxel grid of points (an acceleration structure)
  const double voxel_width = midRadius * 2.0;
  Grid<Eigen::Vector3d> grid(min_bound, max_bound, voxel_width);
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    const Eigen::Vector3d pos = cloud.ends[i];
    grid.insert(grid.index(pos), pos);
  }
  const int min_num_points = 6;

  // 2. initialise one trunk candidate for each occupied voxel
  initialiseTrunks(trunks, cloud, min_bound, voxel_width);
  std::cout << "initialised to " << trunks.size() << " trunks" << std::endl;

  // 3. iterate every candidate several times
  best_trunks_ = trunks;
  for (auto &trunk : best_trunks_)
  {
    trunk.active = false;
  }
  const int num_iterations = 5;
  for (int it = 0; it < num_iterations; it++)
  {
    double above_count = 0;
    double active_count = 0;
    for (int trunk_id = 0; trunk_id < static_cast<int>(trunks.size()); trunk_id++)
    {
      auto &trunk = trunks[trunk_id];
      if (!trunk.active)
      {
        continue;
      }
      // get overlapping points to this trunk
      std::vector<Eigen::Vector3d> points = trunk.getOverlappingPoints(grid, spacing);
      if (points.size() < min_num_points)  // not enough data to use
      {
        trunk.active = false;
        continue;
      }

      // improve the estimation of the trunk's pose and size
      trunk.updateDirection(points);
      trunk.updateCentre(points);
      trunk.updateRadius(points);
      trunk.updateScore(points);

      if (trunk.score > best_trunks_[trunk_id].score)  // got worse, so analyse the best result now
      {
        best_trunks_[trunk_id] = trunk;
      }
      if (trunk.last_score > 0.0 && trunk.score + 3.0 * (trunk.score - trunk.last_score) < minimum_score)
      {
        trunk.active = false;
      }

      bool leaning_too_much = false;
      leaning_too_much = std::abs(best_trunks_[trunk_id].dir[2]) < 0.85;
      if (trunk.length < 4.0 * midRadius || leaning_too_much)  // not enough data to use
      {
        trunk.active = false;
      }

      if (trunk.active)
      {
        active_count++;
      }
      if (trunk.score > minimum_score)
      {
        above_count++;
      }
    }
    std::cout << "iteration " << it << " / " << num_iterations << " " << trunks.size() << " trunks, " << above_count
              << " valid" << std::endl;
    std::cout << active_count << " active" << std::endl;
    if (verbose)
    {
      drawTrunks(trunks, offset);
    }
  }
  trunks.clear();
  // keep only trunks that are above the minimum score, and not leaning too much
  for (auto &trunk : best_trunks_)
  {
    bool leaning_too_much = std::abs(trunk.dir[2]) < 0.9;
    if (trunk.active && trunk.score >= minimum_score && !leaning_too_much)
    {
      trunks.push_back(trunk);
    }
  }
  best_trunks_ = trunks;
  removeOverlappingTrunks(best_trunks_);
  trunks.clear();
  
  for (auto &trunk : best_trunks_)
  {
    if (trunk.dir[2] < 0.0) // orient trunks upwards
    {
      trunk.dir *= -1.0;
    }
    if (trunk.active) // keep only the active trunks
    {
      trunks.push_back(trunk);
    }
  }
  if (verbose)
  {
    drawTrunks(trunks, offset);
    std::cout << "num non-overlapping trunks: " << trunks.size() << std::endl;
  }

  // what if there is a clear cylindrical shape in the points, but rays pass through its centre
  // then it cannot be a trunk. We deal with that situation here
  if (remove_permeable_trunks)
  {
    removePermeableTrunks(verbose, cloud, offset, trunks, min_bound, max_bound);
  }

  setTrunkGroundHeights(cloud, trunks, min_bound, max_bound);
  // find only the lowest trunk to the ground in any near-vertical chain of candidates
  std::vector<int> lowest_trunk_ids = findLowestTrunks(trunks);

  saveDebugTrunks("trunks_verbose.ply", verbose, lowest_trunk_ids, trunks, offset);

  best_trunks_.clear();
  for (auto &id : lowest_trunk_ids)
  {
    best_trunks_.push_back(trunks[id]);
  }
}

/// remove trunk candidates with rays that pass right through them
void Trunks::removePermeableTrunks(bool verbose, const Cloud &cloud, const Eigen::Vector3d &offset, std::vector<Trunk> &trunks, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound)
{
  // first we make a 2D grid to store which horizontal cells (pixels) have rays passing through them
  RayIndexGrid2D grid2D;
  grid2D.init(min_bound, max_bound, 2.0);
  for (auto &trunk : trunks)
  {
    if (trunk.active)
    {
      grid2D.pixel(trunk.centre).filled = true;
    }
  }
  // set the grid's occupancy from the rays
  grid2D.fillRays(cloud);

  std::vector<Eigen::Vector3d> closest_approach_points, pass_through_points;
  int num_removed = 0;
  // now check how occupied the pixels are that overlap each trunk
  for (auto &trunk : trunks)
  {
    if (!trunk.active)
    {
      continue;
    }

    Eigen::Vector3d base = trunk.centre - trunk.length * 0.5 * trunk.dir;
    auto &ray_ids = grid2D.pixel(trunk.centre).ray_ids;
    double mean_rad = 0.0;
    double mean_num = 0.0;
    std::vector<Eigen::Vector3d> nearest_points;
    for (size_t i = 0; i < ray_ids.size(); i++)
    {
      // check whether ray passes through trunk...
      const Eigen::Vector3d start = cloud.starts[ray_ids[i]];
      const Eigen::Vector3d end = cloud.ends[ray_ids[i]];

      const Eigen::Vector3d line_dir = end - start;
      const Eigen::Vector3d across = line_dir.cross(trunk.dir);
      const Eigen::Vector3d side = across.cross(trunk.dir);

      const double d = std::max(0.0, std::min((base - start).dot(side) / line_dir.dot(side), 1.0));
      const Eigen::Vector3d closest_point = start + line_dir * d;

      const double d2 =
        std::max(0.0, std::min((closest_point - base).dot(trunk.dir), trunk.length + 2.0 * trunk.radius));
      const Eigen::Vector3d closest_point2 = base + trunk.dir * d2;

      const double dist = (closest_point - closest_point2).norm();
      if (dist < trunk.radius)
      {
        mean_rad += dist;
        mean_num++;
        nearest_points.push_back(closest_point);
      }
    }
    if (mean_num > 1)
    {
      mean_rad /= mean_num;
      if (mean_rad < 0.7 * trunk.radius)  // too many pass through the trunk
      {
        trunk.active = false;
        num_removed++;
        pass_through_points.insert(pass_through_points.begin(), nearest_points.begin(), nearest_points.end());
      }
      else
      {
        closest_approach_points.insert(closest_approach_points.begin(), nearest_points.begin(), nearest_points.end());
      }
    }
  }
  if (verbose)
  {
    // visualise the trunks and the removed ones, showing the closest point of passing rays
    std::cout << "num trunks removed: " << num_removed << std::endl;
    drawTrunks(trunks, offset, &closest_approach_points, &pass_through_points);
  }

  best_trunks_ = trunks;
  trunks.clear();
  for (auto &trunk : best_trunks_)
  {
    if (trunk.active)
    {
      trunks.push_back(trunk);
    }
  }
}

/// set the ground heights for each trunk
void Trunks::setTrunkGroundHeights(const Cloud &cloud, std::vector<Trunk> &trunks, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound)
{
  // now get ground height for each trunk:
  const double pixel_width = 2.0;
  const int dimx = static_cast<int>(std::ceil((max_bound[0] - min_bound[0]) / pixel_width));
  const int dimy = static_cast<int>(std::ceil((max_bound[1] - min_bound[1]) / pixel_width));
  Eigen::ArrayXXd ground_heights(dimx, dimy);
  ground_heights.setZero();
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    const Eigen::Vector3i index = ((cloud.ends[i] - min_bound) / pixel_width).cast<int>();
    double &h = ground_heights(index[0], index[1]);
    if (h == 0.0 || cloud.ends[i][2] < h)
    {
      h = cloud.ends[i][2];
    }
  }
  // set the trunk ground height, and orient the trunks to not point downwards
  for (auto &trunk : trunks)
  {
    const Eigen::Vector3i index = ((trunk.centre - min_bound) / pixel_width).cast<int>();
    trunk.ground_height = ground_heights(index[0], index[1]) - 1.0;  // TODO: fix!
  }
}

/// a forest nearest path search to find only the lowest trunks of any connected chain
std::vector<int> Trunks::findLowestTrunks(const std::vector<Trunk> &trunks) const
{
  // Next h
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;
  std::vector<int> lowest_trunk_ids;
  // get the lowest points and fill in closest_node
  for (size_t i = 0; i < trunks.size(); i++)
  {
    size_t j;
    for (j = 0; j < trunks.size(); j++)
    {
      Eigen::Vector3d dif = trunks[j].centre - trunks[i].centre;
      dif[2] = 0.0;
      const double x2 = dif.squaredNorm();
      const double height = std::max(0.001, trunks[j].centre[2] - trunks[j].ground_height);

      const double z = x2 / (2.0 * height);
      const Eigen::Vector3d basei = trunks[i].centre - trunks[i].dir * 0.5 * trunks[i].length;
      const Eigen::Vector3d basej = trunks[j].centre - trunks[j].dir * 0.5 * trunks[j].length;
      if (basei[2] > basej[2] + z)
      {
        break;
      }
    }
    if (j == trunks.size())
    {
      closest_node.push(QueueNode(sqr(trunks[i].centre[2]), static_cast<int>(i)));
      lowest_trunk_ids.push_back(static_cast<int>(i));
    }
  }
  std::cout << "number of ground trunks: " << closest_node.size() << " so " << trunks.size() - closest_node.size()
            << " trunks have been removed" << std::endl;
  return lowest_trunk_ids;
}

/// render trunk points to disk:
void Trunks::saveDebugTrunks(const std::string &filename, bool verbose, const std::vector<int> &lowest_trunk_ids, const std::vector<Trunk> &trunks, const Eigen::Vector3d &offset) const
{
  if (verbose)
  {
    std::vector<Eigen::Vector3d> cloud_points;
    std::vector<double> times;
    std::vector<RGBA> colours;
    RGBA colour;
    colour.alpha = 255;
    for (auto &id : lowest_trunk_ids)
    {
      const Trunk &trunk = trunks[id];
      colour.red = uint8_t(rand() % 255);
      colour.green = uint8_t(rand() % 255);
      colour.blue = uint8_t(rand() % 255);

      const Eigen::Vector3d side1 = trunk.dir.cross(Eigen::Vector3d(1, 2, 3)).normalized();
      const Eigen::Vector3d side2 = side1.cross(trunk.dir);
      const double rad = trunk.radius;
      for (double z = -0.5; z < 0.5; z += 0.1)
      {
        for (double ang = 0.0; ang < 2.0 * kPi; ang += 0.3)
        {
          const Eigen::Vector3d pos =
            trunk.centre + trunk.dir * trunk.length * z + side1 * std::sin(ang) * rad + side2 * std::cos(ang) * rad;
          cloud_points.push_back(pos + offset);
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }
    writePlyPointCloud(filename, cloud_points, times, colours);
  }  
}

bool Trunks::save(const std::string &filename, const Eigen::Vector3d &offset) const
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << "# tree trunks file:" << std::endl;
  ofs << "x,y,z,radius" << std::endl;
  for (auto &trunk : best_trunks_)
  {
    if (!trunk.active)
    {
      continue;
    }
    Eigen::Vector3d base = trunk.centre - trunk.dir * trunk.length * 0.5;
    ofs << base[0]+offset[0] << ", " << base[1]+offset[1] << ", " << base[2]+offset[2] << ", " << trunk.radius << std::endl;
  }
  return true;
}

std::vector<std::pair<Eigen::Vector3d, double>> Trunks::load(const std::string &filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for reading." << std::endl;
    return std::vector<std::pair<Eigen::Vector3d, double>>();
  }
  std::vector<std::pair<Eigen::Vector3d, double>> trunks;
  for (std::string line; std::getline(ifs, line);)
  {
    Eigen::Vector3d base;
    double radius;
    if (line.length() == 0 || line[0] == '#')
    {
      continue;
    }
    const int num_commas = static_cast<int>(std::count(line.begin(), line.end(), ','));
    if (num_commas == 3)  // just the base
    {
      std::istringstream ss(line);
      for (int i = 0; i < 4; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i < 3)
        {
          base[i] = std::stod(token.c_str());
        }
        else
        {
          radius = std::stod(token.c_str());
        }
      }
      trunks.push_back(std::pair<Eigen::Vector3d, double>(base, radius));
    }
    else
    {
      std::cerr << "bad input, there should be 4 fields per line: x, y, z, radius." << std::endl;
      return std::vector<std::pair<Eigen::Vector3d, double>>();
    }
  }
  return trunks;
}


}  // namespace ray
