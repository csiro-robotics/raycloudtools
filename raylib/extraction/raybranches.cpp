// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raybranches.h"
#include "../raydebugdraw.h"
#include "../raygrid.h"
#include "../raycuboid.h"
#include <map>

namespace ray
{
#define LEAN
namespace
{
struct Branch
{
  Branch() : centre(0,0,0), radius(0), score(0), length(0), dir(0,0,0) {}
  Eigen::Vector3d centre; 
  double radius;
  double score;
  double length; 
  Eigen::Vector3d dir;
};

struct Accumulator
{
  Accumulator(): weight(0), x(0), abs_x(0), y(0,0), xy(0,0), x2(0), radius(0), radius2(0) {}

  double weight;
  double x;
  double abs_x;
  Eigen::Vector2d y;
  Eigen::Vector2d xy;
  double x2;
  double radius;
  double radius2;
};

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

void drawBranches(const std::vector<Branch> &branches)
{
  std::vector<Eigen::Vector3d> starts(branches.size()), ends(branches.size());
  std::vector<double> radii(branches.size());
  for (size_t i = 0; i<branches.size(); i++)
  {
    starts[i] = branches[i].centre - branches[i].dir*branches[i].length*0.5;
    ends[i] = branches[i].centre + branches[i].dir*branches[i].length*0.5;
    radii[i] = branches[i].radius;
  }
  DebugDraw::instance()->drawCylinders(starts, ends, radii, 1);
}

static const double branch_height_to_width = 4.0; // height extent relative to real diameter of branch
static const double boundary_radius_scale = 2.0; // how much farther out is the expected boundary compared to real branch radius? Larger requires more space to declare it a branch

void getOverlap(const Grid<Eigen::Vector3d> &grid, const Branch &branch, std::vector<Eigen::Vector3d> &points)
{
  Eigen::Vector3d base = branch.centre - 0.5*branch.length*branch.dir;
  Eigen::Vector3d top = branch.centre + 0.5*branch.length*branch.dir;
  double outer_radius = branch.radius * boundary_radius_scale;
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
        for (auto &pos: cell.data)
        {
          double h = pos[2] - branch.centre[2];
          if (std::abs(h) > branch.length*0.5)
          {
            continue;
          }
          double dist2 = (branch.centre - pos).squaredNorm();
          if (dist2 <= outer_radius*outer_radius)
          {
            points.push_back(pos);
          }
        }
      }
    }
  }
}
}

Bush::Bush(const Cloud &cloud, double midRadius, bool verbose)
{
  double spacing = cloud.estimatePointSpacing();
  
  // Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
  const double minimum_score = 50.0; 
  
  if (verbose)
  {
    std::cout << "av radius: " << midRadius << ", estimated point spacing: " << spacing << ", minimum score: " << minimum_score << std::endl;
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
  const int min_num_points = 6;

  // 2. initialise one branch candidate for each occupied voxel
  std::vector<Branch> branches;
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
  auto addBranch = [min_num_points, &branches]
    (double voxel_width, const Eigen::Vector3d &offset, const Eigen::Vector3i &index, int count)
  {
    Eigen::Vector3d centre = (index.cast<double>() + Eigen::Vector3d(0.5,0.5,0.5))*voxel_width + offset;
    if (count < 2)
      return;
    Branch branch;
    double diameter = voxel_width / std::sqrt(2.0);
    branch.centre = centre;
    branch.radius = diameter / 2.0;
    branch.length = diameter * branch_height_to_width;
    branch.score = 0;
    branch.dir = Eigen::Vector3d(0,0,1);
    branches.push_back(branch);    
  };
  for (int j = 0; j<4; j++)
  {
    voxels[j]->forEach(addBranch);
  }
  if (verbose)
  {
    drawBranches(branches);
  }  
  // 3. iterate every candidate several times
  const int num_iterations = 5;
  for (int it = 0; it<num_iterations; it++)
  {
    for (int branch_id = 0; branch_id < (int)branches.size(); branch_id++)
    {
      auto &branch = branches[branch_id];
      // get overlapping points to this branch
      std::vector<Eigen::Vector3d> points;
      getOverlap(grid, branch, points);
      if (points.size() < min_num_points) // not enough data to use
      {
        branches[branch_id] = branches.back(); // so remove the branch
        branches.pop_back();
        branch_id--;
        continue;
      }

      // get an initial direction vector estimate if there isn't one already
      if (branch.dir[2] == 1.0)
      {
        Eigen::Vector3d centroid = mean(points);
        Eigen::Matrix3d scatter;
        scatter.setZero();
        for (auto &point: points)
          scatter += (point - centroid) * (point - centroid).transpose();
        scatter /= (double)points.size();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());   
        branch.dir = eigen_solver.eigenvectors().col(2);   
  //      branch.radius = std::sqrt(eigen_solver.eigenvalues()[0]) + std::sqrt(eigen_solver.eigenvalues()[1]);
        branch.centre = centroid;
        continue;
      }

      // 1. estimate branch direction, and distance up the branch
      {
        Eigen::Vector3d ax1 = Eigen::Vector3d(1,2,3).cross(branch.dir).normalized();
        Eigen::Vector3d ax2 = ax1.cross(branch.dir);
        Accumulator sum;
        for (size_t i = 0; i<points.size(); i++)
        {
          Eigen::Vector3d to_point = points[i] - branch.centre;
          Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));

          double dist = offset.norm();
          double w = 1.0 - dist/(branch.radius * boundary_radius_scale); // lateral fade off
          // remove radius. If radius_removal_factor=0 then half-sided trees will have estimated branch centred on that edge
          //                If radius_removal_factor=1 then v thin branches may accidentally get a radius and it won't shrink down
          const double radius_removal_factor = 0.5;
          offset -= offset * radius_removal_factor * branch.radius / offset.norm(); 

          double h = to_point.dot(branch.dir);
          // lean, shift and change radius
          sum.x += h*w;
          sum.y += offset*w;
          sum.xy += h*offset*w;
          sum.x2 += h*h*w;
          sum.abs_x += std::abs(h)*w;
          sum.weight += w;      
        }
        double n = sum.weight;
        branch.centre += branch.dir*(sum.x / n);

        // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
        Eigen::Vector2d sXY = sum.xy - sum.x*sum.y/n;
        double sXX = sum.x2 - sum.x*sum.x/n;
        if (std::abs(sXX) > 1e-10)
          sXY /= sXX;

        branch.dir = (branch.dir + ax1*sXY[0] + ax2*sXY[1]).normalized();
//        branch.length = 4.0*(sum.abs_x/n);
      }

      // 2. get new centre of cylinder
      {
        Eigen::Vector3d ax1 = Eigen::Vector3d(1,2,3).cross(branch.dir).normalized();
        Eigen::Vector3d ax2 = ax1.cross(branch.dir);

        std::vector<Eigen::Vector3d> ps(points.size());
        Eigen::Vector3d mean_p(0,0,0);
        for (size_t i = 0; i<points.size(); i++)
        {
          Eigen::Vector3d to_point = points[i] - branch.centre;
          Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));
          Eigen::Vector2d xy = offset/branch.radius;
          double l2 = xy.squaredNorm();
          Eigen::Vector3d point(xy[0], xy[1], 0.5*l2); // a paraboloid that has gradient 1 at 1
          ps[i] = point;
          mean_p += point;
        }
        mean_p /= (double)points.size();      
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
        double A = (plane.xz*plane.y2 - plane.yz*plane.xy) / (plane.x2*plane.y2 - plane.xy*plane.xy);
        double B = (plane.yz - A * plane.xy) / plane.y2;
        Eigen::Vector2d shift(A,B);
        double shift2 = shift.squaredNorm();
        if (shift2 > 1.0) // don't shift more than one radius each iteration
          shift /= std::sqrt(shift2);

        branch.centre += (ax1*shift[0] + ax2*shift[1]) * branch.radius;   
      }

      // 3. estimate radius and score
      {
        Accumulator sum;
        branch.score = 0;
        std::vector<double> scores(points.size());
        for (size_t i = 0; i<points.size(); i++)
        {
          Eigen::Vector3d to_point = points[i] - branch.centre;
          double dist_sqr = (to_point - branch.dir*branch.dir.dot(to_point)).squaredNorm();
          sum.radius += std::sqrt(dist_sqr);
          sum.radius2 += dist_sqr;
        }
        double n = (double)points.size();
        branch.radius = sum.radius / n;
        branch.length = 2.0*branch.radius * branch_height_to_width;
        double num_points = (double)points.size() - 4.0; // (double)min_num_points;
        double variance = (sum.radius2/n - sqr(sum.radius / n)) * n/num_points; // end part gives sample variance
        double density = num_points * sqr(spacing) / (2.0 * kPi * branch.radius * branch.length);
        branch.score = std::sqrt(density / variance);
    
        if (it == num_iterations-1 && branch.score < minimum_score) // then remove the branch
        {
          branches[branch_id] = branches.back(); 
          branches.pop_back();
          branch_id--;
          continue;        
        }
      }
      if (branch.radius > 0.5*branch.length || branch.length < midRadius) // not enough data to use
      {
        branches[branch_id] = branches.back(); // so remove the branch
        branches.pop_back();
        branch_id--;
        continue;
      }
    }
    if (verbose)
    {
      drawBranches(branches);
      std::cout << "num branches: " << branches.size() << std::endl;
    }  
  }
  if (verbose)
  {
    drawBranches(branches);
  } 


  // Next a forest nearest path search

}

bool Bush::save(const std::string &filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }  
  ofs << "# Tree base location list: x, y, z, radius" << std::endl;
 /* for (auto &branch: branch_bases)
  {
    Eigen::Vector3d base = branch.centre - branch.dir*branch.length*0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << branch.radius << std::endl;
  }*/
  return true;
}

std::vector<std::pair<Eigen::Vector3d, double> > Bush::load(const std::string &filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for reading." << std::endl;
    return std::vector<std::pair<Eigen::Vector3d, double> >();
  }  
  std::vector<std::pair<Eigen::Vector3d, double> > branches;
  while (!ifs.eof())
  {
    Eigen::Vector3d base;
    double radius;
    std::string line;
    std::getline(ifs, line);
    if (line.length() == 0 || line[0] == '#')
      continue;
    int num_commas = (int)std::count(line.begin(), line.end(), ',');
    if (num_commas == 3) // just the base
    {
      std::istringstream ss(line);
      for (int i = 0; i<4; i++)
      {
        std::string token;
        std::getline(ss, token, ',');
        if (i<3)
          base[i] = std::stod(token.c_str());
        else
          radius = std::stod(token.c_str());
      }
      branches.push_back(std::pair<Eigen::Vector3d, double>(base, radius));
    }
    else
    {
      std::cerr << "bad input, there should be 4 fields per line: x, y, z, radius." << std::endl;
      return std::vector<std::pair<Eigen::Vector3d, double> >();
    }
  }
  return branches;  
}


} // namespace ray
