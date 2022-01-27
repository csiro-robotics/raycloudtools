// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrunks.h"
#include "../raydebugdraw.h"
#include "../raygrid.h"
#include "../raycuboid.h"
#include "../rayforestgen.h"
#include <map>

namespace ray
{
#define LEAN
namespace
{
struct Ray
{
  Ray(const Eigen::Vector3d &start, const Eigen::Vector3d &pos) : start(start), pos(pos) {}
  Eigen::Vector3d start, pos;
};

struct Cell
{
  std::vector<Ray> rays;
  Eigen::Vector2d minBound;
  double height;
};

struct Accumulator
{
  Accumulator(): weight(0), x(0), abs_x(0), y(0,0), xy(0,0), x2(0), radius(0), radius2(0), z(0,0), xz(0,0) {}

  double weight;
  double x;
  double abs_x;
  Eigen::Vector2d y;
  Eigen::Vector2d xy;
  double x2;
  double radius;
  double radius2;
  Eigen::Vector2d z;
  Eigen::Vector2d xz;
};

Eigen::Vector2d vector2d(const Eigen::Vector3d &v)
{
  return Eigen::Vector2d(v[0], v[1]);
}
Eigen::Vector3d vector3d(const Eigen::Vector2d &v, double z = 0)
{
  return Eigen::Vector3d(v[0], v[1], z);
}

struct IntegerVoxels
{
  IntegerVoxels(double width, const Eigen::Vector3d offset) : voxel_width(width), offset(offset) {}

  inline Eigen::Vector3i getIndex(const Eigen::Vector3d &pos)
  {
    Eigen::Vector3d ind = (pos - offset) / voxel_width;
    return Eigen::Vector3i(int(std::floor(ind[0])), int(std::floor(ind[1])), int(std::floor(ind[2])));
  }
  inline void increment(const Eigen::Vector3d &pos)
  {
    increment(getIndex(pos));
  }
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
  inline int get(const Eigen::Vector3i &index)
  {
    auto it = voxel_map.find(index);
    if (it == voxel_map.end())
      return 0;
    else
      return it->second;
  }
  void forEach(std::function<void(double width, const Eigen::Vector3d &offset, const Eigen::Vector3i &index, int count)> func)
  {
    for (auto &voxel: voxel_map)
    {
      func(voxel_width, offset, voxel.first, voxel.second);
    }
  }

  std::map<Eigen::Vector3i, int, Vector3iLess> voxel_map;
  double voxel_width;
  Eigen::Vector3d offset;
};

static const double trunk_height_to_width = 3.0; // height extent relative to real radius of trunk
static const double boundary_radius_scale = 2.0; // how much farther out is the expected boundary compared to real trunk radius? Larger requires more space to declare it a trunk

void getOverlap(const Grid<Eigen::Vector3d> &grid, const Trunk &trunk, std::vector<Eigen::Vector3d> &points)
{
  Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
  Eigen::Vector3d base = trunk.centre - 0.5*trunk.length*lean;
  Eigen::Vector3d top = trunk.centre + 0.5*trunk.length*lean;
  double outer_radius = trunk.radius * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, 0);
  Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
//  std::cout << "radius: " << trunk.radius << ", min: " << cuboid.min_bound_.transpose() << ", max: " << cuboid.max_bound_.transpose() << std::endl;
//  std::cout << "box min: " << grid.box_min.transpose() << std::endl;
//  std::cout << "mins: " << mins.transpose() << ", maxs: " << maxs.transpose() << std::endl;
  // just in case...
  mins = maxVector(mins, Eigen::Vector3i(0,0,0));
  Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1,1,1);
  maxs = minVector(maxs, min_dims);

  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0]<=maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1]<=maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2]<=maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        for (auto &pos: cell.data)
        {
          double h = pos[2] - trunk.centre[2];
          if (std::abs(h) > trunk.length*0.5)
          {
            continue;
          }
          Eigen::Vector3d centre = trunk.centre + lean * h;
          double dist2 = (centre - pos).squaredNorm();
          if (dist2 <= outer_radius*outer_radius)
          {
            points.push_back(pos);
          }
        }
      }
    }
  }
}


void getOverlap(const Grid<Trunk> &grid, const Trunk &trunk, std::vector<Trunk*> &overlapping_trunks)
{
  Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
  const double length_scale = 2.0; // make trunks longer when searching for overlapping trunks
  Eigen::Vector3d base = trunk.centre - 0.5*trunk.length*length_scale*lean;
  Eigen::Vector3d top = trunk.centre + 0.5*trunk.length*length_scale*lean;
  double outer_radius = trunk.radius * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, 0);
  Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();

  mins = maxVector(mins, Eigen::Vector3i(0,0,0));
  Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1,1,1);
  maxs = minVector(maxs, min_dims);

  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0]<=maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1]<=maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2]<=maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        for (auto &other_trunk: cell.data)
        {
          if (&other_trunk == &trunk)
            continue;
          // cylinder-cylinder overlap
          double lowest = std::max(trunk.centre[2]-0.5*trunk.length*length_scale, other_trunk.centre[2]-0.5*other_trunk.length*length_scale);
          double highest = std::min(trunk.centre[2]+0.5*trunk.length*length_scale, other_trunk.centre[2]+0.5*other_trunk.length*length_scale);
          if (highest <= lowest)
            continue;
          #if 1 // better
          Eigen::Vector3d oth = other_trunk.centre + vector3d(other_trunk.lean, 1)*(trunk.centre[2]-other_trunk.centre[2]) - trunk.centre;
          Eigen::Vector3d oth_dir = vector3d(other_trunk.lean, 1) - vector3d(trunk.lean, 1);
          oth_dir[2] = 0;
          double d = std::max(lowest-trunk.centre[2], std::min(oth.dot(oth_dir)/oth_dir.dot(oth_dir), highest-trunk.centre[2]));
          Eigen::Vector3d intersection = oth + oth_dir*d;
          double dist = intersection.norm();
          #else
          double mid = (lowest + highest)/2.0;
          Eigen::Vector3d p1 = trunk.centre + vector3d(trunk.lean, 1)*(mid - trunk.centre[2]);
          Eigen::Vector3d p2 = other_trunk.centre + vector3d(other_trunk.lean, 1)*(mid - other_trunk.centre[2]);
          Eigen::Vector3d p = (p2 - p1);
          double dist = p.norm();
          #endif
          if (dist > 1.25*(trunk.radius + other_trunk.radius))
            continue;
          overlapping_trunks.push_back((Trunk *)&other_trunk);
        }
      }
    }
  }
}
}

Wood::Wood(const Cloud &cloud, double midRadius, bool verbose)
{
  double spacing = cloud.estimatePointSpacing();
  
  // Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
  // trunk thickness affects how strictly it adheres to a cylinder.
  // lower trunk_thickness (stricter) will require a lower minimum_score to find the same number of trees.
  // at the point where minimum score is 0, it is invariant to the number of points
  #define SIGMA_SCORE
  #if defined SIGMA_SCORE
    const double minimum_score = 50.0; 
  #else
    const double minimum_score = 0.3/sqr(spacing);
  #endif
 // const double minimum_score = 2000.0;
  
  if (verbose)
  {
    std::cout << "estimated point spacig: " << spacing << ", minimum score: " << minimum_score << std::endl;
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 0);
  }

  Eigen::Vector3d min_bound = cloud.calcMinBound();
  Eigen::Vector3d max_bound = cloud.calcMaxBound();
  std::cout << "cloud from: " << min_bound.transpose() << " to: " << max_bound.transpose() << std::endl;

  // 1. voxel grid of points (an acceleration structure)
  const double voxel_width = midRadius * 2.0;
  Grid<Eigen::Vector3d> grid(min_bound, max_bound, voxel_width);
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3d pos = cloud.ends[i];
    grid.insert(grid.index(pos), pos);
  }
  const int min_num_points = 10;

  // 2. initialise one trunk candidate for each occupied voxel
  std::vector<Trunk> trunks;
  #if 1
  Eigen::Vector3d half_voxel(0.5*voxel_width, 0.5*voxel_width, 0.5*voxel_width);
  // normal size, with offset grid
  IntegerVoxels voxels1(voxel_width, min_bound), voxels2(voxel_width, min_bound + half_voxel);
  // double size with offset grid
  IntegerVoxels voxels3(2.0*voxel_width, min_bound), voxels4(2.0*voxel_width, min_bound + 2.0*half_voxel);
  IntegerVoxels *voxels[4] = {&voxels1, &voxels2, &voxels3, &voxels4};
  Eigen::Vector3i ind2;
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    for (int j = 0; j<4; j++)
    {
      voxels[j]->increment(cloud.ends[i]);
    }
  }
  auto addTrunk = [min_num_points, &trunks]
    (double width, const Eigen::Vector3d &offset, const Eigen::Vector3i &index, int count)
  {
    Eigen::Vector3d centre = (index.cast<double>() + Eigen::Vector3d(0.5,0.5,0.5))*width + offset;
    if (count < min_num_points)
      return;
    Trunk trunk;
    double diameter = width / std::sqrt(2.0);
    trunk.centre = centre;
    trunk.radius = diameter / 2.0;
    trunk.length = diameter * trunk_height_to_width;
    trunk.score = trunk.weight = trunk.combined_score = 0;
    trunk.lean.setZero();
    trunks.push_back(trunk);    
  };
  for (int j = 0; j<4; j++)
  {
    voxels[j]->forEach(addTrunk);
  }
  #else
  grid.walkCells([&](const Grid<Eigen::Vector3d> &grid, const Grid<Eigen::Vector3d>::Cell &cell)
  {
    if (cell.data.size() < min_num_points) // not enough data to attempt
      return;
    Trunk trunk;
    trunk.centre = (cell.index.cast<double>() + Eigen::Vector3d(0.5,0.5,0.5)) * voxel_width + grid.box_min;
    trunk.radius = midRadius;
    trunk.length = 2.0 * midRadius * trunk_height_to_width;
    trunk.score = trunk.weight = trunk.combined_score = 0;
    trunk.lean.setZero();
    trunks.push_back(trunk);
  });
  #endif
  if (verbose)
  {
    DebugDraw::instance()->drawTrunks(trunks);
  }  

  // 3. iterate every candidate several times
  const int num_iterations = 5;
  for (int it = 0; it<num_iterations; it++)
  {
    for (int trunk_id = 0; trunk_id < (int)trunks.size(); trunk_id++)
    {
      auto &trunk = trunks[trunk_id];
      // get overlapping points to this trunk
      std::vector<Eigen::Vector3d> points;
      getOverlap(grid, trunk, points);
      if (points.size() < min_num_points) // not enough data to use
      {
        trunks[trunk_id] = trunks.back(); // so remove the trunk
        trunks.pop_back();
        trunk_id--;
        continue;
      }

      // weight these points
      Eigen::Vector3d lean(trunk.lean[0], trunk.lean[1], 1);
      std::vector<double> weights(points.size());

      Accumulator sum;
      trunk.score = 0;
      std::vector<double> scores(points.size());
      std::vector<Eigen::Vector3d> ps(points.size());
      struct Acc
      {
        Acc(){ x2 = y2 = xy = xz = yz = 0; }
        double x2, y2, xy, xz, yz;
      };

      Eigen::Vector3d mean_p(0,0,0);
      for (size_t i = 0; i<points.size(); i++)
      {
        double h = points[i][2] - trunk.centre[2];
        Eigen::Vector2d offset = vector2d(points[i] - (trunk.centre + lean*h));
        Eigen::Vector2d xy = offset/trunk.radius;
        double l2 = xy.squaredNorm();
        Eigen::Vector3d point(xy[0], xy[1], 0.5*l2); // a paraboloid that has gradient 1 at 1
        ps[i] = point;
        mean_p += point;

        double dist = offset.norm();
        double w = 1.0 - dist/(trunk.radius * boundary_radius_scale); // lateral fade off
        weights[i] = w;
        // remove radius. If radius_removal_factor=0 then half-sided trees will have estimated trunk centred on that edge
        //                If radius_removal_factor=1 then v thin trunks may accidentally get a radius and it won't shrink down
        const double radius_removal_factor = 0.5;
        offset -= offset * radius_removal_factor * trunk.radius / offset.norm(); 

        // lean, shift and change radius
        sum.x += h*w;
        sum.y += offset*w;
        sum.xy += h*offset*w;
        sum.x2 += h*h*w;
        sum.radius += dist;
        sum.radius2 += dist*dist;
        sum.abs_x += std::abs(h)*w;
        sum.weight += w;      

#if !defined SIGMA_SCORE
        const double trunk_thickness = 0.025; 
        // hard coding for now. Representing the expected error from circular in metres for real trees
        double score_centre = 1.0 - trunk.radius/trunk_thickness;
        double score_radius = 1.0;
        double score_2radius = 1.0 - trunk.radius/trunk_thickness;
        double score_3radius = 1.0 - trunk.radius/trunk_thickness;
        double weight = 0.0;
        if (dist < trunk.radius)
          weight = score_centre + (score_radius - score_centre)*dist/trunk.radius;
        else if (dist < 2.0*trunk.radius)
          weight = score_radius + (score_2radius - score_radius)*(dist - trunk.radius)/trunk.radius;
        else
          weight = score_2radius + (score_3radius - score_2radius)*(dist - 2.0*trunk.radius)/trunk.radius;

        if (it == num_iterations-1)
        {
          scores[i] = weight;
        }
        trunk.score += weight / (2.0 * kPi * trunk.radius * trunk.length);
#endif
      }
      double N = (double)points.size();
      mean_p /= N;
#if defined SIGMA_SCORE
      double num_points = (double)points.size() - 4.0; // (double)min_num_points;
      double variance = (sum.radius2/N - sqr(sum.radius / N)) * N/num_points; // end part gives sample variance
      double density = num_points * sqr(spacing) / (2.0 * kPi * trunk.radius * trunk.length);
      trunk.score = std::sqrt(density / variance);
#endif    
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
      double n = sum.weight;

      if (it == num_iterations-1 && trunk.score < minimum_score) // then remove the trunk
      {
        trunks[trunk_id] = trunks.back(); 
        trunks.pop_back();
        trunk_id--;
        continue;        
      }

      // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
      Eigen::Vector2d sXY = sum.xy - sum.x*sum.y/n;
      double sXX = sum.x2 - sum.x*sum.x/n;

      trunk.lean += sXY / sXX;
      double l = trunk.lean.norm();
      const double max_lean = 0.5;//25;
      if (l > max_lean)
        trunk.lean *= max_lean/l;          
      
      double A = (plane.xz*plane.y2 - plane.yz*plane.xy) / (plane.x2*plane.y2 - plane.xy*plane.xy);
      double B = (plane.yz - A * plane.xy) / plane.y2;

      Eigen::Vector2d shift(A,B);
      double height = mean_p[2] + (shift[0]-mean_p[0])*A + (shift[1]-mean_p[1])*B;
      double paraboloid_height = 0.5*(A*A + B*B);
      double new_radius = std::sqrt(2.0*(height - paraboloid_height)) * trunk.radius;

      new_radius = std::min(new_radius, trunk.radius * 2.0); // don't grow by more than a factor of 2 per iteration.
      double shift2 = shift.squaredNorm();
      if (shift2 > 1.0) // don't shift more than one radius each iteration
        shift /= std::sqrt(shift2);

      trunk.centre += vector3d(shift * trunk.radius);   

      double radius_scale = new_radius / trunk.radius;
 //     double radius_scale = (sum.radius/(double)points.size()) / trunk.radius;

      trunk.centre[2] += sum.x / n;
      double length_scale = (sum.abs_x/n) / (trunk.length * 0.25);
      trunk.radius *= radius_scale;
      trunk.length *= length_scale;
      if (trunk.radius > 0.5*trunk.length || trunk.length < midRadius) // not enough data to use
      {
        trunks[trunk_id] = trunks.back(); // so remove the trunk
        trunks.pop_back();
        trunk_id--;
        continue;
      }
    }
    if (verbose)
    {
      DebugDraw::instance()->drawTrunks(trunks);
    }  
  }
  if (verbose)
  {
    DebugDraw::instance()->drawTrunks(trunks);
  } 

  // now I need to connect all the trunks into tree shapes
  Grid<Trunk> trunk_grid(min_bound, max_bound, voxel_width);
  for (auto &trunk: trunks)
  {
    trunk.next_down = NULL;
    trunk_grid.insert(trunk_grid.index(trunk.centre), trunk);
  }
  for (size_t i = 0; i<trunks.size(); i++)
  {
    std::vector<Trunk*> overlaps;
    getOverlap(trunk_grid, trunks[i], overlaps);
    for (auto &overlap: overlaps)
    {
      // no, not quite, we need an insertion sort type thing...
      if (trunks[i].centre[2] > overlap->centre[2])
      {
        if (trunks[i].next_down == NULL || trunks[i].next_down->centre[2] < overlap->centre[2])
          trunks[i].next_down = overlap;
      }
      else
      {
        if (overlap->next_down == NULL || overlap->next_down->centre[2] < trunks[i].centre[2])
          overlap->next_down = &trunks[i];
      }
    }
  }
  std::vector<Trunk *> trunk_pointers;
  for (auto &trunk: trunks) // TODO: we can also walk this and generate a tree of branches
  {
    Trunk *tr = &trunk;
    while (tr->next_down != NULL)
      tr = tr->next_down;
    if (tr->combined_score == 0)
      trunk_pointers.push_back(tr);
    tr->combined_score += trunk.score;
  }
  const double consensus_scale = 2.0; // requires this many minimum scores to pass. Higher requires longer tree shapes
  for (auto &trunk: trunk_pointers)
  {
    if (trunk->combined_score > minimum_score * consensus_scale)
    {
      trunk_bases.push_back(*trunk);
    }
  }
  double mean_radius = 0;
  double total_volume = 0;
  const double volume_gain = 35.0; // this is the factor for the artificial trees in raycreate
  for (auto &t: trunk_bases)
  {
    mean_radius += t.radius;
    total_volume += volume_gain * t.radius * t.radius * t.radius;
  }
  mean_radius /= (double)trunk_bases.size();
  std::cout << "number of trees found: " << trunk_bases.size() << " with mean diameter: " << 2.0*mean_radius << std::endl;
  std::cout << "approximate volume of wood: " << total_volume << " cubic metres" << std::endl;
  const double wood_density = 500; // for Red gum, see: https://www.engineeringtoolbox.com/wood-density-d_40.html
  std::cout << "estimated mass of trees (at 500kg/m^3): " << total_volume * wood_density / 1000.0 << " tonnes" << std::endl;
  if (verbose)
  {
    DebugDraw::instance()->drawTrunks(trunk_bases);
  } 
}

bool Wood::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# trunks file:" << std::endl;
  ofs << "x,y,z,radius" << std::endl;
  for (auto &trunk: trunk_bases)
  {
    Eigen::Vector3d base = trunk.centre - vector3d(trunk.lean, 1)*trunk.length*0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << trunk.radius << std::endl;
  }
  return true;
}

std::vector<std::pair<Eigen::Vector3d, double> > Wood::load(const std::string &filename)
{
  ForestStructure forest;
  std::vector<std::pair<Eigen::Vector3d, double> > trunks;
  forest.load(filename);
  for (auto &tree: forest.trees)
  {
    trunks.push_back(std::pair<Eigen::Vector3d, double>(tree.segments()[0].tip, tree.segments()[0].radius));
  }
  return trunks;  
}


} // namespace ray
