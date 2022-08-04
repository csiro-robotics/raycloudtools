// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrunks.h"
#include <nabo/nabo.h>
#include <queue>
#include "../raycuboid.h"
#include "../raydebugdraw.h"
#include "../raygrid.h"
#include "../rayply.h"

namespace ray
{
namespace
{
const double inf = 1e10;

void drawTrunks(const std::vector<Trunk> &trunks, const std::vector<Eigen::Vector3d> *extra_points1 = NULL,
                const std::vector<Eigen::Vector3d> *extra_points2 = NULL)
{
#if USE_RVIZ
  std::vector<Eigen::Vector3d> starts, ends;
  std::vector<double> radii;
  std::vector<Eigen::Vector3d> colours;
  for (size_t i = 0; i < trunks.size(); i++)
  {
    if (!trunks[i].active)
      continue;
    if (!(trunks[i].radius == trunks[i].radius))
      std::cout << "nan trunk radius at " << i << std::endl;
    if (trunks[i].radius <= 0.0)
      std::cout << "bad trunk radius at " << i << std::endl;
    if (!(trunks[i].centre == trunks[i].centre))
      std::cout << "nan trunk centre at " << i << std::endl;
    starts.push_back(trunks[i].centre - trunks[i].dir * trunks[i].length * 0.5);
    ends.push_back(trunks[i].centre + trunks[i].dir * trunks[i].length * 0.5);
    radii.push_back(trunks[i].radius);
    double shade = std::min(trunks[i].score / (2.0 * minimum_score), 1.0);
    Eigen::Vector3d col(0, 0, 0);
    col[0] = col[1] = shade;
    col[2] = shade > 0.5 ? 1.0 : 0.0;
    colours.push_back(col);
  }
  //  DebugDraw::instance()->drawCylinders(starts, ends, radii, 1, colours);
  DebugDraw::instance()->drawLines(starts, ends, colours);
#else
  std::vector<Eigen::Vector3d> cloud_points;
  std::vector<double> times;
  std::vector<RGBA> colours;
  RGBA colour;
  colour.red = colour.green = colour.blue = 0;
  colour.alpha = 255;

  if (extra_points1)
  {
    for (auto &point : *extra_points1)
    {
      cloud_points.push_back(point);
      times.push_back(0.0);
      colours.push_back(colour);
    }
  }
  colour.red = colour.green = 255;
  if (extra_points2)
  {
    for (auto &point : *extra_points2)
    {
      cloud_points.push_back(point);
      times.push_back(0.0);
      colours.push_back(colour);
    }
  }
  for (auto &trunk : trunks)
  {
    if (!trunk.active)
      continue;
    colour.red = colour.green = colour.blue = (uint8_t)std::min(255.0 * trunk.score / (2.0 * minimum_score), 255.0);
    if (trunk.score < minimum_score)
      colour.green = colour.blue = 0;

    Eigen::Vector3d side1 = trunk.dir.cross(Eigen::Vector3d(1, 2, 3)).normalized();
    Eigen::Vector3d side2 = side1.cross(trunk.dir);
    double rad = trunk.radius;
    for (double z = -0.5; z < 0.5; z += 0.15)
    {
      for (double ang = 0.0; ang < 2.0 * kPi; ang += 0.6)
      {
        Eigen::Vector3d pos =
          trunk.centre + trunk.dir * trunk.length * z + side1 * std::sin(ang) * rad + side2 * std::cos(ang) * rad;
        cloud_points.push_back(pos);
        times.push_back(0.0);
        colours.push_back(colour);
      }
    }
  }
  writePlyPointCloud("trunks_verbose.ply", cloud_points, times, colours);
#endif
}

/// Used in a priority queue for identifying close neighbouring trunks
struct QueueNode
{
  QueueNode() {}
  QueueNode(double score, int index)
    : score(score)
    , id(index)
  {}

  double score;
  int id;
};

//#define MINIMISE_ANGLE // works quite well in flowing along trunks, but sometimes causes multi-trunk problem, where
//radius was too small.
class QueueNodeComparator
{
public:
  bool operator()(const QueueNode &p1, const QueueNode &p2) { return p1.score > p2.score; }
};
}  // namespace


/// Initial estimate of the trunks
void initialiseTrunks(std::vector<Trunk> &trunks, const Cloud &cloud, const Eigen::Vector3d &min_bound,
                      double voxel_width)
{
  // we generate a set of candidate trunks at multiple resolutions
  Eigen::Vector3d half_voxel(0.5 * voxel_width, 0.5 * voxel_width, 0.5 * voxel_width);
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
      continue;
    for (int j = 0; j < 4; j++)
    {
      voxels[j]->increment(cloud.ends[i]);
    }
  }
  auto addTrunk = [&trunks](const struct IntegerVoxels &voxels, double voxel_width, const Eigen::Vector3d &offset,
                            const Eigen::Vector3i &index, int count) {
    Eigen::Vector3d centre = (index.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width + offset;
    if (count < 4) // not enough points to add a candidate trunk
      return;
    // now check that there are points above and below the voxel in question. To make sure that there could be
    // a tall enough trunk here in the data
    const int minpoints = 4;
    bool above1 = voxels.get(index + Eigen::Vector3i(0, 0, 1)) >= minpoints;
    bool above2 = voxels.get(index + Eigen::Vector3i(0, 0, 2)) >= minpoints;
    bool below1 = voxels.get(index - Eigen::Vector3i(0, 0, 1)) >= minpoints;
    bool below2 = voxels.get(index - Eigen::Vector3i(0, 0, 2)) >= minpoints;
    bool tall_enough = (above1 && below1) || (above1 && above2) || (below1 && below2);
    if (!tall_enough)
      return;
    // ok this is good enough for a candidate, so add it in
    Trunk trunk;
    double diameter = voxel_width / std::sqrt(2.0);
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
void removeOverlappingTrunks(std::vector<Trunk> &best_trunks)
{
  // Brute force approach at the moment!
  for (size_t i = 0; i < best_trunks.size(); i++)
  {
    if (!(i % 1000))
      std::cout << i << " / " << best_trunks.size() << std::endl;
    Trunk &trunk = best_trunks[i];
    if (!trunk.active)
      continue;
    
    // ideally we would get the volume of overlap of two cylinders, however
    // since this is complicated, we approximate by taking a multiple points in the 
    // first cylinder and counting how many overlap the second cylinder
    Eigen::Vector3d ax1 = Eigen::Vector3d(1, 2, 3).cross(trunk.dir).normalized();
    Eigen::Vector3d ax2 = trunk.dir.cross(ax1);
    const int num = 5;
    double s = 0.8;
    // 5 points in a cross 
    double xs[num] = { 0, s, 0, -s, 0 };
    double ys[num] = { 0, 0, s, 0, -s };
    // x 5 heights along the cylinder = 25 points
    double zs[num] = { -0.5 * s, -0.25 * s, 0, 0.25 * s, 0.5 * s };

    // for coarse intersection
    Eigen::Vector3d base = trunk.centre - 0.5 * trunk.length * trunk.dir;
    Eigen::Vector3d top = trunk.centre + 0.5 * trunk.length * trunk.dir;
    Eigen::Vector3d rad(trunk.radius, trunk.radius, trunk.radius);
    Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

    for (size_t j = 0; j < best_trunks.size(); j++) // brute force
    {
      if (i == j)
        continue;
      Trunk &cylinder = best_trunks[j];
      if (!cylinder.active)
        continue;
      // first, bounding box intersection test
      Eigen::Vector3d base2 = cylinder.centre - 0.5 * cylinder.length * cylinder.dir;
      Eigen::Vector3d top2 = cylinder.centre + 0.5 * cylinder.length * cylinder.dir;
      Eigen::Vector3d rad2(cylinder.radius, cylinder.radius, cylinder.radius);
      Cuboid cuboid2(minVector(base2, top2) - rad2, maxVector(base2, top2) + rad2);
      if (!cuboid.overlaps(cuboid2))  // broadphase exclusion
        continue;

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
          double d = pos.dot(cylinder.dir);
          if (d > cylinder.length * 0.5 || d < -cylinder.length * 0.5)
            continue;
          pos -= cylinder.dir * d;
          if (pos.squaredNorm() < sqr(cylinder.radius))
          {
            num_inside++;
          }
        }
      }

      double inside_ratio = (double)num_inside / (double)(num * num);
      if (inside_ratio > 0.4)
      {
        double vol_trunk = sqr(trunk.radius) * trunk.length;
        double vol_cylinder = sqr(cylinder.radius) * cylinder.length;
        // remove the smaller one
        if (vol_trunk < vol_cylinder)
          trunk.active = false;
        else
          cylinder.active = false;
        break;
      }
    }
  }
}

// trunk identification in ray cloud
Trunks::Trunks(const Cloud &cloud, double midRadius, bool verbose, bool exclude_passing_rays)
{
  // spacing helps us choose a grid resolution that is appropriate to the cloud
  double spacing = cloud.estimatePointSpacing();  
  if (verbose)
  {
    std::cout << "av radius: " << midRadius << ", estimated point spacing: " << spacing
              << ", minimum score: " << minimum_score << std::endl;
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 1);
  }
  Eigen::Vector3d min_bound = cloud.calcMinBound();
  Eigen::Vector3d max_bound = cloud.calcMaxBound();
  std::cout << "cloud from: " << min_bound.transpose() << " to: " << max_bound.transpose() << std::endl;

  std::vector<Trunk> trunks;

  // 1. voxel grid of points (an acceleration structure)
  const double voxel_width = midRadius * 2.0;
  Grid<Eigen::Vector3d> grid(min_bound, max_bound, voxel_width);
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3d pos = cloud.ends[i];
    grid.insert(grid.index(pos), pos);
  }
  const int min_num_points = 6;

  // 2. initialise one trunk candidate for each occupied voxel
  initialiseTrunks(trunks, cloud, min_bound, voxel_width);
  std::cout << "initialised to " << trunks.size() << " trunks" << std::endl;

  RayGrid2D grid2D;
  grid2D.init(min_bound, max_bound, 2.0);

  // 3. iterate every candidate several times
  best_trunks = trunks;
  for (auto &trunk : best_trunks) trunk.active = false;
  const int num_iterations = 5;
  for (int it = 0; it < num_iterations; it++)
  {
    double above_count = 0;
    double active_count = 0;
    for (int trunk_id = 0; trunk_id < (int)trunks.size(); trunk_id++)
    {
      auto &trunk = trunks[trunk_id];
      if (!trunk.active)
        continue;
      // get overlapping points to this trunk
      std::vector<Eigen::Vector3d> points;
      trunk.getOverlap(grid, points, spacing);
      if (points.size() < min_num_points)  // not enough data to use
      {
        trunk.active = false;
        continue;
      }

      // improve the estimation of the trunk's pose and size
      trunk.updateDirection(points);
      trunk.updateCentre(points);
      trunk.updateRadiusAndScore(points);

      if (trunk.score > best_trunks[trunk_id].score)  // got worse, so analyse the best result now
      {
        best_trunks[trunk_id] = trunk;
      }
      if (trunk.last_score > 0.0 && trunk.score + 3.0 * (trunk.score - trunk.last_score) < minimum_score)
      {
        trunk.active = false;
      }

      bool leaning_too_much = false;
      leaning_too_much = std::abs(best_trunks[trunk_id].dir[2]) < 0.85;
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
      drawTrunks(trunks);
    }
  }
  trunks.clear();
  // keep only trunks that are above the minimum score, and not leaning too much
  for (auto &trunk : best_trunks)
  {
    bool leaning_too_much = false;
    leaning_too_much = std::abs(trunk.dir[2]) < 0.9;
    if (trunk.active && trunk.score >= minimum_score && !leaning_too_much)  
    {
      trunks.push_back(trunk);
    }
  }
  best_trunks = trunks;
  removeOverlappingTrunks(best_trunks);
  trunks.clear();
  // keep only the active trunks
  for (auto &trunk : best_trunks)
  {
    if (trunk.active)
    {
      trunks.push_back(trunk);
    }
  }
  if (verbose)
  {
    drawTrunks(trunks);
    std::cout << "num non-overlapping trunks: " << trunks.size() << std::endl;
  }

  // what if there is a clear cylindrical shape in the points, but rays pass through its centre
  // then it cannot be a trunk. We deal with that situation here
  if (exclude_passing_rays)
  {  
    for (auto &trunk : trunks)
    {
      if (trunk.active)
      {
        grid2D.pixel(trunk.centre).filled = true;
      }
    }
    // set occupancy from ray cloud
    grid2D.fillRays(cloud);

    std::vector<Eigen::Vector3d> extra_points1, extra_points2;
    int num_removed = 0;
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
      std::vector<Eigen::Vector3d> poins;
      for (size_t i = 0; i < ray_ids.size(); i++)
      {
        // check whether ray passes through trunk...
        Eigen::Vector3d start = cloud.starts[ray_ids[i]];
        Eigen::Vector3d end = cloud.ends[ray_ids[i]];

        Eigen::Vector3d line_dir = end - start;
        Eigen::Vector3d across = line_dir.cross(trunk.dir);
        Eigen::Vector3d side = across.cross(trunk.dir);

        double d = std::max(0.0, std::min((base - start).dot(side) / line_dir.dot(side), 1.0));
        Eigen::Vector3d closest_point = start + line_dir * d;

        double d2 = std::max(0.0, std::min((closest_point - base).dot(trunk.dir), trunk.length + 2.0 * trunk.radius));
        Eigen::Vector3d closest_point2 = base + trunk.dir * d2;

        double dist = (closest_point - closest_point2).norm();
        if (dist < trunk.radius)
        {
          mean_rad += dist;
          mean_num++;
          poins.push_back(closest_point);
        }
      }
      if (mean_num > 1)
      {
        mean_rad /= mean_num;
        if (mean_rad < 0.7 * trunk.radius) // too many pass through the trunk
        {
          trunk.active = false;
          num_removed++;
          extra_points2.insert(extra_points2.begin(), poins.begin(), poins.end());
        }
        else
        {
          extra_points1.insert(extra_points1.begin(), poins.begin(), poins.end());
        }
      }
    }
    if (verbose)
    {
      // visualise the trunks and the removed ones, showing the closest point of passing rays
      std::cout << "num trunks removed: " << num_removed << std::endl;
      drawTrunks(trunks, &extra_points1, &extra_points2);
    }

    best_trunks = trunks;
    trunks.clear();
    for (auto &trunk : best_trunks)
    {
      if (trunk.active)
        trunks.push_back(trunk);
    }
  }

  // get lowest trunks:
  // now get ground height for each trunk:
  double pixel_width = 2.0;
  int dimx = (int)std::ceil((max_bound[0] - min_bound[0]) / pixel_width);
  int dimy = (int)std::ceil((max_bound[1] - min_bound[1]) / pixel_width);
  Eigen::ArrayXXd ground_heights(dimx, dimy);
  ground_heights.setZero();
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3i index = ((cloud.ends[i] - min_bound) / pixel_width).cast<int>();
    double &h = ground_heights(index[0], index[1]);
    if (h == 0.0 || cloud.ends[i][2] < h)
      h = cloud.ends[i][2];
  }
  for (auto &trunk : trunks)
  {
    Eigen::Vector3i index = ((trunk.centre - min_bound) / pixel_width).cast<int>();
    trunk.ground_height = ground_heights(index[0], index[1]) - 1.0;  // TODO: fix!
    if (trunk.dir[2] < 0.0)
      trunk.dir *= -1.0;
  }

  // Next a forest nearest path search
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
      double x2 = dif.squaredNorm();
      double height = std::max(0.001, trunks[j].centre[2] - trunks[j].ground_height);

      double z = x2 / (2.0 * height);
      Eigen::Vector3d basei = trunks[i].centre - trunks[i].dir * 0.5 * trunks[i].length;
      Eigen::Vector3d basej = trunks[j].centre - trunks[j].dir * 0.5 * trunks[j].length;
      if (basei[2] > basej[2] + z)
      {
        break;
      }
    }
    if (j == trunks.size())
    {
      if (trunks[i].dir[2] < 0.0)
      {
        trunks[i].dir *= -1;
      }
      closest_node.push(QueueNode(sqr(trunks[i].centre[2]), (int)i));
      lowest_trunk_ids.push_back((int)i);
    }
  }
  std::cout << "number of ground trunks: " << closest_node.size() << " so " << trunks.size() - closest_node.size()
            << " trunks have been removed" << std::endl;

  // render these trunk points to disk:
  if (verbose)
  {
    std::vector<Eigen::Vector3d> cloud_points;
    std::vector<double> times;
    std::vector<RGBA> colours;
    RGBA colour;
    colour.alpha = 255;
    for (auto &id : lowest_trunk_ids)
    {
      Trunk &trunk = trunks[id];
      colour.red = uint8_t(rand() % 255);
      colour.green = uint8_t(rand() % 255);
      colour.blue = uint8_t(rand() % 255);

      Eigen::Vector3d side1 = trunk.dir.cross(Eigen::Vector3d(1, 2, 3)).normalized();
      Eigen::Vector3d side2 = side1.cross(trunk.dir);
      double rad = trunk.radius;
      for (double z = -0.5; z < 0.5; z += 0.1)
      {
        for (double ang = 0.0; ang < 2.0 * kPi; ang += 0.3)
        {
          Eigen::Vector3d pos =
            trunk.centre + trunk.dir * trunk.length * z + side1 * std::sin(ang) * rad + side2 * std::cos(ang) * rad;
          cloud_points.push_back(pos);
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }
    writePlyPointCloud("trunks_verbose.ply", cloud_points, times, colours);
  }

  best_trunks.clear();
  for (auto &id : lowest_trunk_ids)
  {
    best_trunks.push_back(trunks[id]);
  }
}

bool Trunks::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << "# tree trunks file:" << std::endl;
  ofs << "x,y,z,radius" << std::endl;
  for (auto &trunk : best_trunks)
  {
    if (!trunk.active)
      continue;
    Eigen::Vector3d base = trunk.centre - trunk.dir * trunk.length * 0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << trunk.radius << std::endl;
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
  while (!ifs.eof())
  {
    Eigen::Vector3d base;
    double radius;
    std::string line;
    std::getline(ifs, line);
    if (line.length() == 0 || line[0] == '#')
      continue;
    int num_commas = (int)std::count(line.begin(), line.end(), ',');
    if (num_commas == 3)  // just the base
    {
      std::istringstream ss(line);
      for (int i = 0; i < 4; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i < 3)
          base[i] = std::stod(token.c_str());
        else
          radius = std::stod(token.c_str());
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
