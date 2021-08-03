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
#include <nabo/nabo.h>
#include <queue>

namespace ray
{
#define LEAN
namespace
{
const double inf = 1e10;
struct Branch
{
  Branch() : centre(0,0,0), radius(0), score(0), length(0), dir(0,0,0), parent(-1), tree_score(inf), distance_to_ground(inf), active(true), visited(false) {}
  Eigen::Vector3d centre; 
  double radius;
  double score;
  double length; 
  Eigen::Vector3d dir;
  int parent;
  double tree_score;
  double distance_to_ground;
  bool active;
  bool visited;
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

// Tuning: minimum_score defines how sparse your tree feature can be, compared to the decimation spacing
static const double minimum_score = 40.0; 
static const double branch_height_to_width = 4.0; // height extent relative to real diameter of branch
static const double boundary_radius_scale = 2.0; // how much farther out is the expected boundary compared to real branch radius? Larger requires more space to declare it a branch

void drawBranches(const std::vector<Branch> &branches)
{
  std::vector<Eigen::Vector3d> starts(branches.size()), ends(branches.size());
  std::vector<double> radii(branches.size());
  std::vector<Eigen::Vector4d> colours(branches.size());
  for (size_t i = 0; i<branches.size(); i++)
  {
    if (!branches[i].active)
      continue;
    starts.push_back(branches[i].centre - branches[i].dir*branches[i].length*0.5);
    ends.push_back(branches[i].centre + branches[i].dir*branches[i].length*0.5);
    radii.push_back(branches[i].radius);
    double shade = std::min(branches[i].score / (2.0 * minimum_score), 1.0);
    Eigen::Vector4d col(0,0,0,0.5);
    col[0] = col[1] = shade;
    col[2] = shade > 0.5 ? 1.0 : 0.0;
    colours.push_back(col);
  }
  DebugDraw::instance()->drawCylinders(starts, ends, radii, 1, colours);
}

struct QueueNode
{
  QueueNode(){}
  QueueNode(double distance_to_ground, double score, int index) : distance_to_ground(distance_to_ground), score(score), id(index) {}

  double distance_to_ground;
  double score;
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


void getOverlap(const Grid<Eigen::Vector3d> &grid, const Branch &branch, std::vector<Eigen::Vector3d> &points, double spacing)
{
  Eigen::Vector3d base = branch.centre - 0.5*branch.length*branch.dir;
  Eigen::Vector3d top = branch.centre + 0.5*branch.length*branch.dir;
  double outer_radius = (branch.radius + spacing) * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, outer_radius);
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
          Eigen::Vector3d p = pos - branch.centre;
          double h = p.dot(branch.dir);
          if (std::abs(h) > branch.length*0.5)
          {
            continue;
          }
          p -= branch.dir*h;
          double dist2 = p.squaredNorm();
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
  
  if (verbose)
  {
    std::cout << "av radius: " << midRadius << ", estimated point spacing: " << spacing << ", minimum score: " << minimum_score << std::endl;
    DebugDraw::instance()->drawCloud(cloud.ends, 0.5, 1);
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
//    drawBranches(branches);
  }  
  // 3. iterate every candidate several times
  std::vector<Branch> best_branches = branches;
  for (auto &branch: best_branches)
    branch.active = false;
  const int num_iterations = 5;
  for (int it = 0; it<num_iterations; it++)
  {
    std::cout << "iteration " << it << " / " << num_iterations << " " << branches.size() << " branches" << std::endl;
    for (int branch_id = 0; branch_id < (int)branches.size(); branch_id++)
    {
      auto &branch = branches[branch_id];
      if (!branch.active)
        continue;
      // get overlapping points to this branch
      std::vector<Eigen::Vector3d> points;
      getOverlap(grid, branch, points, spacing);
      if (points.size() < min_num_points) // not enough data to use
      {
        branch.active = false;
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
          // const double dist = offset.norm();
          double w = 1.0;//1.0 - dist/(branch.radius * boundary_radius_scale); // lateral fade off
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
     //   branch.centre += branch.dir*(sum.x / n); // in theory it moves towards a better spot, but in practice it gets rid of diversity of positions

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
        if (branch.score > best_branches[branch_id].score) // got worse, so analyse the best result now
          best_branches[branch_id] = branch;      
      }
      if (branch.radius > 0.5*branch.length || branch.length < midRadius) // not enough data to use
        branch.active = false;
    }
  }
  if (verbose)
  {
//    drawBranches(best_branches);
  }   
  for (int branch_id = 0; branch_id<(int)best_branches.size(); branch_id++)
  {
    if (!best_branches[branch_id].active || best_branches[branch_id].score < minimum_score) // then remove the branch
    {
      best_branches[branch_id] = best_branches.back(); 
      best_branches.pop_back();
      branch_id--;
      continue;        
    }    
  }
  if (verbose)
  {
 //   drawBranches(best_branches);
  } 
  std::cout << "num valid branches: " << best_branches.size() << std::endl;
  // Next, clean up the set of branches by removing overlapping ones
  // (brute force for now)
  for (size_t i = 0; i<best_branches.size(); i++)
  {
    if (!(i%1000))
      std::cout << i << " / " << best_branches.size() << std::endl;
    Branch &branch = best_branches[i];
    if (!branch.active)
      continue;
    Eigen::Vector3d ax1 = Eigen::Vector3d(1,2,3).cross(branch.dir).normalized();
    Eigen::Vector3d ax2 = branch.dir.cross(ax1);
    const int num = 5;
    double s = 0.8;
    double xs[num] = {0, s, 0, -s, 0};
    double ys[num] = {0, 0, s,  0, -s};
    double zs[num] = {-0.5*s, -0.25*s, 0, 0.25*s, 0.5*s};

    // for coarse intersection
    Eigen::Vector3d base = branch.centre - 0.5*branch.length*branch.dir;
    Eigen::Vector3d top = branch.centre + 0.5*branch.length*branch.dir;
    Eigen::Vector3d rad(branch.radius, branch.radius, branch.radius);
    Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);    

    for (size_t j = 0; j<best_branches.size(); j++)
    {
      if (i==j)
        continue;
      Branch &cylinder = best_branches[j];
      if (!cylinder.active)
        continue;
      Eigen::Vector3d base2 = cylinder.centre - 0.5*cylinder.length*cylinder.dir;
      Eigen::Vector3d top2 = cylinder.centre + 0.5*cylinder.length*cylinder.dir;
      Eigen::Vector3d rad2(cylinder.radius, cylinder.radius, cylinder.radius);
      Cuboid cuboid2(minVector(base2, top2) - rad2, maxVector(base2, top2) + rad2);  
      if (!cuboid.overlaps(cuboid2)) // broadphase exclusion
        continue;

      int num_inside = 0;
      for (int k = 0; k<num; k++)
      {
        for (int l = 0; l<num; l++)
        {
          Eigen::Vector3d pos = branch.centre + branch.dir*zs[k]*branch.length + (ax1*xs[l] + ax2*ys[l])*branch.radius;
          pos -= cylinder.centre;
          // is pos inside cylinder j?
          double d = pos.dot(cylinder.dir);
          if (d > cylinder.length*0.5 || d < -cylinder.length*0.5)
            continue;
          pos -= cylinder.dir*d;
          if (pos.squaredNorm() < sqr(cylinder.radius))
            num_inside++;
        }
      }
      double inside_ratio = (double)num_inside / (double)(num*num);
      if (inside_ratio > 0.4)
      {
        double vol_branch = sqr(branch.radius)*branch.length;
        double vol_cylinder = sqr(cylinder.radius)*cylinder.length;
        // remove the smaller one
        if (vol_branch < vol_cylinder)
        {
          best_branches[i] = best_branches.back();
          best_branches.pop_back();
          i--;
        }
        else
        {
          if (j > i)
          {
            best_branches[j] = best_branches.back();
            best_branches.pop_back();
          }
          else 
            cylinder.active = false;
        }
        break;
      }
    }
  }
  branches.clear();
  for (auto &branch: best_branches)
    if (branch.active)
      branches.push_back(branch);
  std::cout << "num non-overlapping branches: " << branches.size() << std::endl;
  
  if (verbose)
  {
    drawBranches(branches);
  } 

  // Next a forest nearest path search

	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;
  
  // get the lowest points and fill in closest_node
  for (size_t i = 0; i<branches.size(); i++)
  {
    size_t j;
    for (j = 0; j<branches.size(); j++)
    {
      Eigen::Vector3d dif = branches[j].centre - branches[i].centre;
      dif[2] = 0.0;
      double x2 = dif.squaredNorm();
      double height = branches[j].centre[2] - min_bound[2]; // TODO: use ground map for better value here

      double z = x2/(2.0*height);
      if (branches[i].centre[2] > z)
        break;
    }
    if (j==branches.size())
    {
      closest_node.push(QueueNode(branches[i].centre[2], sqr(branches[i].centre[2]), (int)i));
    }
  }
  std::cout << "number of ground branches: " << closest_node.size() << std::endl;

  // 1. get nearest neighbours
  const int search_size = 20;
  Eigen::MatrixXd points_p(3, branches.size());
  for (unsigned int i = 0; i < branches.size(); i++) 
    points_p.col(i) = branches[i].centre;
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, branches.size());
  dists2.resize(search_size, branches.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0);
  
  // 2b. climb up from lowest branches...
	while(!closest_node.empty())
  {
		QueueNode node = closest_node.top(); closest_node.pop();
		if(!branches[node.id].visited)
    {
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        double dist = std::sqrt(dists2(i, node.id));
        double new_dist = node.distance_to_ground + dist;///branches[node.id].radius;
        double new_score = 0;
        #if defined MINIMISE_SQUARE_DISTANCE
        dist *= dist;
        #endif
        #if defined MINIMISE_ANGLE
        Eigen::Vector3d dif = (branches[child].centre - branches[node.id].centre).normalized();
        Eigen::Vector3d dir = branches[node.id].dir;
        Eigen::Vector3d dir2 = branches[child].dir;
        if (dir2.dot(dir) < 0.0)
          dir2 = -dir2;
        dir = (dir + dir2).normalized();
        
        const double power = 2.0;
        dist /= std::pow(std::max(0.001, dif.dot(dir)), power);
        #endif
     //   dist /= branches[node.id].radius;
        #if defined MINIMISE_SQUARE_DISTANCE || defined MINIMISE_ANGLE
        new_score = node.score + dist;
        if (new_score < branches[child].tree_score)
        #else
        if (new_dist < points[child].distance_to_ground)
        #endif
        {
					branches[child].tree_score = new_score;
					branches[child].distance_to_ground = new_dist;
          branches[child].parent = node.id;
					closest_node.push(QueueNode(branches[child].distance_to_ground, branches[child].tree_score, child));
				}
			}
		  branches[node.id].visited = true;
		}
	}
  // now we need to render the structure as a tree... I'll use lines
  std::vector<Eigen::Vector3d> starts, ends, colours;
  for (auto &branch: branches)
  {
    if (branch.parent == -1)
      continue;
    starts.push_back(branch.centre);
    ends.push_back(branches[branch.parent].centre);
    Eigen::Vector3d col;
    col[0] = std::fmod(branch.tree_score,       1.0);
    col[1] = std::fmod(branch.tree_score/10.0,  1.0);
    col[2] = std::fmod(branch.tree_score/100.0, 1.0);
    colours.push_back(col);
  }
  DebugDraw::instance()->drawLines(starts, ends, colours);
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
