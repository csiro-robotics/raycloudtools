// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raybranches.h"
#include "../raydebugdraw.h"
#include "../raygrid.h"
#include "../raycuboid.h"
#include "../rayply.h"
#include <nabo/nabo.h>
#include <queue>

namespace ray
{
namespace
{
const double inf = 1e10;

void drawBranches(const std::vector<Branch> &branches, const std::vector<Eigen::Vector3d> *extra_points1 = NULL, const std::vector<Eigen::Vector3d> *extra_points2 = NULL)
{
  #if USE_RVIZ
  std::vector<Eigen::Vector3d> starts, ends;
  std::vector<double> radii;
  std::vector<Eigen::Vector3d> colours;
  for (size_t i = 0; i<branches.size(); i++)
  {
    if (!branches[i].active)
      continue;
    if (!(branches[i].radius == branches[i].radius))
      std::cout << "nan branch radius at " << i << std::endl;
    if (branches[i].radius <= 0.0)
      std::cout << "bad branch radius at " << i << std::endl;
    if (!(branches[i].centre == branches[i].centre))
      std::cout << "nan branch centre at " << i << std::endl;
    starts.push_back(branches[i].centre - branches[i].dir*branches[i].length*0.5);
    ends.push_back(branches[i].centre + branches[i].dir*branches[i].length*0.5);
    radii.push_back(branches[i].radius);
    double shade = std::min(branches[i].score / (2.0 * minimum_score), 1.0);
    Eigen::Vector3d col(0,0,0);
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
      for (auto &point: *extra_points1)
      {
        cloud_points.push_back(point);
        times.push_back(0.0);
        colours.push_back(colour);
      }
    }
    colour.red = colour.green = 255;
    if (extra_points2)
    {
      for (auto &point: *extra_points2)
      {
        cloud_points.push_back(point);
        times.push_back(0.0);
        colours.push_back(colour);
      }
    }    
    for (auto &branch: branches)
    {
      if (!branch.active)
        continue;
  //    colour.red = uint8_t(rand()%255);
  //    colour.green = uint8_t(rand()%255);
  //    colour.blue = uint8_t(rand()%255);
      colour.red = colour.green = colour.blue = (uint8_t)std::min(255.0*branch.score / (2.0*minimum_score), 255.0);
      if (branch.score < minimum_score)
        colour.green = colour.blue = 0;

      Eigen::Vector3d side1 = branch.dir.cross(Eigen::Vector3d(1,2,3)).normalized();
      Eigen::Vector3d side2 = side1.cross(branch.dir);
      double rad = branch.radius;
      for (double z = -0.5; z<0.5; z+=0.15)
      {
        for (double ang = 0.0; ang<2.0*kPi; ang += 0.6)
        {
          Eigen::Vector3d pos = branch.centre + branch.dir*branch.length*z + side1*std::sin(ang)*rad + side2*std::cos(ang)*rad;
          cloud_points.push_back(pos);
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }
    writePlyPointCloud("branches_verbose.ply", cloud_points, times, colours);
  #endif
}

struct QueueNode
{
  QueueNode(){}
  QueueNode(double score, int index) : score(score), id(index) {}

  double score;
  int id;
};

//#define MINIMISE_ANGLE // works quite well in flowing along branches, but sometimes causes multi-branch problem, where radius was too small. 
class QueueNodeComparator 
{ 
public: 
    bool operator() (const QueueNode &p1, const QueueNode &p2) 
    { 
      return p1.score > p2.score; 
    } 
}; 
}

void initialiseBranches(std::vector<Branch> &branches, const Cloud &cloud, const Eigen::Vector3d &min_bound, double voxel_width)
{
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
  auto addBranch = [&branches]
    (const struct IntegerVoxels &voxels, double voxel_width, const Eigen::Vector3d &offset, const Eigen::Vector3i &index, int count)
  {
    Eigen::Vector3d centre = (index.cast<double>() + Eigen::Vector3d(0.5,0.5,0.5))*voxel_width + offset;
    if (count < 4)
      return;
    const int minpoints = 4;
    bool above1 = voxels.get(index + Eigen::Vector3i(0,0,1)) >= minpoints;
    bool above2 = voxels.get(index + Eigen::Vector3i(0,0,2)) >= minpoints;
    bool below1 = voxels.get(index - Eigen::Vector3i(0,0,1)) >= minpoints;
    bool below2 = voxels.get(index - Eigen::Vector3i(0,0,2)) >= minpoints;
    bool tall_enough = (above1 && below1) || (above1 && above2) || (below1 && below2);
    if (!tall_enough)
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
}

void removeOverlappingBranches(std::vector<Branch> &best_branches)
{
  // Brute force approach at the moment!
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
      std::vector<Eigen::Vector3d> ps;
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
          branch.active = false;
        else
          cylinder.active = false;
        break;
      }
    }
  }
}

Bush::Bush(const Cloud &cloud, double midRadius, bool verbose, bool trunks_only, bool exclude_passing_rays)
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

  std::vector<Branch> branches;
  
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
  initialiseBranches(branches, cloud, min_bound, voxel_width);
  std::cout << "initialised to " << branches.size() << " branches" << std::endl;

  RayGrid2D grid2D;
  grid2D.init(min_bound, max_bound, 2.0);
 
  // 3. iterate every candidate several times
  best_branches = branches;
  for (auto &branch: best_branches)
    branch.active = false;
  const int num_iterations = 5;
  for (int it = 0; it<num_iterations; it++)
  {
    double above_count = 0;
    double active_count = 0;
    for (int branch_id = 0; branch_id < (int)branches.size(); branch_id++)
    {
      auto &branch = branches[branch_id];
      if (!branch.active)
        continue;
      // get overlapping points to this branch
      std::vector<Eigen::Vector3d> points;
      branch.getOverlap(grid, points, spacing);
      if (points.size() < min_num_points) // not enough data to use
      {
        branch.active = false;
        continue;
      }

      if (branch.dir[2] == 1.0 && !trunks_only)
      {
        branch.estimatePose(points);
        continue;
      }

      branch.updateDirection(points, trunks_only);
      branch.updateCentre(points);
      branch.updateRadiusAndScore(points, spacing, trunks_only);

      if (branch.score > best_branches[branch_id].score) // got worse, so analyse the best result now
        best_branches[branch_id] = branch;      
      if (branch.last_score > 0.0 && branch.score + 3.0*(branch.score - branch.last_score) < minimum_score)
        branch.active = false;

      bool leaning_too_much = false;
      if (trunks_only)
        leaning_too_much = std::abs(best_branches[branch_id].dir[2]) < 0.85;
      if (branch.length < 4.0*midRadius || leaning_too_much) // not enough data to use
        branch.active = false;

      if (branch.active)
        active_count++;
      if (branch.score > minimum_score)
        above_count++;
    }
    std::cout << "iteration " << it << " / " << num_iterations << " " << branches.size() << " branches, " << above_count << " valid" << std::endl;
    std::cout << active_count << " active" << std::endl;
    if (verbose)
    {
      drawBranches(branches); // best_branches);
    } 
  }
  branches.clear();
  for (auto &branch: best_branches)
  {
    bool leaning_too_much = false;
    if (trunks_only)
      leaning_too_much = std::abs(branch.dir[2]) < 0.9;
    if (branch.active && branch.score >= minimum_score && !leaning_too_much) // then remove the branch
      branches.push_back(branch);       
  }  
  best_branches = branches;
  removeOverlappingBranches(best_branches);
  branches.clear();
  for (auto &branch: best_branches)
  {
    if (branch.active) // then remove the branch
      branches.push_back(branch);         
  }    
  if (verbose)
  {
    drawBranches(branches);
    std::cout << "num non-overlapping branches: " << branches.size() << std::endl;
  }   

  if (exclude_passing_rays)
  { // add in internal points due to passing rays...
    for (auto &branch: branches)
    {
      if (branch.active)
        grid2D.pixel(branch.centre).filled = true;
    }
    grid2D.fillRays(cloud); 

    std::vector<Eigen::Vector3d> extra_points1, extra_points2;
    int num_removed = 0;
    for (auto &branch: branches)
    {
      if (!branch.active)
        continue;

      Eigen::Vector3d base = branch.centre - branch.length*0.5*branch.dir;
      auto &ray_ids = grid2D.pixel(branch.centre).ray_ids;
      double mean_rad = 0.0;
      double mean_num = 0.0;
      std::vector<Eigen::Vector3d> poins;
      for (size_t i = 0; i<ray_ids.size(); i++)
      {
        Eigen::Vector3d start = cloud.starts[ray_ids[i]];
        Eigen::Vector3d end = cloud.ends[ray_ids[i]];
        
        Eigen::Vector3d line_dir = end - start;
        Eigen::Vector3d across = line_dir.cross(branch.dir);
        Eigen::Vector3d side = across.cross(branch.dir);

        double d = std::max(0.0, std::min((base - start).dot(side) / line_dir.dot(side), 1.0));
        Eigen::Vector3d closest_point = start + line_dir * d;

        double d2 = std::max(0.0, std::min((closest_point - base).dot(branch.dir), branch.length + 2.0*branch.radius));
        Eigen::Vector3d closest_point2 = base + branch.dir * d2;
        
        double dist = (closest_point - closest_point2).norm();
        if (dist < branch.radius)
        {
          mean_rad += dist;
          mean_num++;
          poins.push_back(closest_point);
        }
      }
      if (mean_num > 1)
      {
        mean_rad /= mean_num;
        if (mean_rad < 0.7*branch.radius)
        {
          branch.active = false;
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
      std::cout << "num trunks removed: " << num_removed << std::endl;
      drawBranches(branches, &extra_points1, &extra_points2);
    }

    best_branches = branches;
    branches.clear();
    for (auto &branch: best_branches)
    {
      if (branch.active) 
        branches.push_back(branch);         
    }    
  }

  // get lowest branches:
  // now get ground height for each branch:
  double pixel_width = 2.0;
  int dimx = (int)std::ceil((max_bound[0] - min_bound[0])/pixel_width);
  int dimy = (int)std::ceil((max_bound[1] - min_bound[1])/pixel_width);
  Eigen::ArrayXXd ground_heights(dimx, dimy);
  ground_heights.setZero();
  for (size_t i = 0; i<cloud.ends.size(); i++)
  {
    if (!cloud.rayBounded(i))
      continue;
    Eigen::Vector3i index = ((cloud.ends[i] - min_bound)/pixel_width).cast<int>();
    double &h = ground_heights(index[0], index[1]);
    if (h == 0.0 || cloud.ends[i][2] < h)
      h = cloud.ends[i][2];
  }
  for (auto &branch: branches)
  {
    Eigen::Vector3i index = ((branch.centre - min_bound)/pixel_width).cast<int>();
    branch.ground_height = ground_heights(index[0], index[1]) - 1.0;  // TODO: fix!
    if (branch.dir[2] < 0.0)
      branch.dir *= -1.0;
  }

  // Next a forest nearest path search
	std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node;
  std::vector<int> lowest_branch_ids;
  // get the lowest points and fill in closest_node
  for (size_t i = 0; i<branches.size(); i++)
  {
    size_t j;
    for (j = 0; j<branches.size(); j++)
    {
      Eigen::Vector3d dif = branches[j].centre - branches[i].centre;
      dif[2] = 0.0;
      double x2 = dif.squaredNorm();
      double height = std::max(0.001, branches[j].centre[2] - branches[j].ground_height); 

      double z = x2/(2.0*height);
      Eigen::Vector3d basei = branches[i].centre - branches[i].dir*0.5*branches[i].length;
      Eigen::Vector3d basej = branches[j].centre - branches[j].dir*0.5*branches[j].length;
      if (basei[2] > basej[2] + z)
      {
//        std::cout << "branch " << i << ": " << basei.transpose() << " is above branch " << j << ": " << basej.transpose() << " with ground height: " << branches[j].ground_height << " so removing" << std::endl;
        break;
      }
    }
    if (j==branches.size())
    {
      if (branches[i].dir[2] < 0.0)
        branches[i].dir *= -1;
      closest_node.push(QueueNode(sqr(branches[i].centre[2]), (int)i));
      lowest_branch_ids.push_back((int)i);
    }
  }
  std::cout << "number of ground trunks: " << closest_node.size() << " so " << branches.size() - closest_node.size() << " trunks have been removed" << std::endl;

  // render these trunk points to disk:
  if (verbose)
  {
    std::vector<Eigen::Vector3d> cloud_points;
    std::vector<double> times;
    std::vector<RGBA> colours;
    RGBA colour;
    colour.alpha = 255;
    for (auto &id: lowest_branch_ids)
    {
      Branch &branch = branches[id];
      colour.red = uint8_t(rand()%255);
      colour.green = uint8_t(rand()%255);
      colour.blue = uint8_t(rand()%255);

      Eigen::Vector3d side1 = branch.dir.cross(Eigen::Vector3d(1,2,3)).normalized();
      Eigen::Vector3d side2 = side1.cross(branch.dir);
      double rad = branch.radius;
      for (double z = -0.5; z<0.5; z+=0.1)
      {
        for (double ang = 0.0; ang<2.0*kPi; ang += 0.3)
        {
          Eigen::Vector3d pos = branch.centre + branch.dir*branch.length*z + side1*std::sin(ang)*rad + side2*std::cos(ang)*rad;
          cloud_points.push_back(pos);
          times.push_back(0.0);
          colours.push_back(colour);
        }
      }
    }
    writePlyPointCloud("trunks_verbose.ply", cloud_points, times, colours);
  }


  if (trunks_only)
  {
    best_branches.clear();
    for (auto &id: lowest_branch_ids)
    {
      best_branches.push_back(branches[id]);
    }
    // now save
    return;
  }
  // 1. get nearest neighbours
  const int search_size = std::min(20, (int)branches.size()-1);
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
      int par = branches[node.id].parent;
      if (par != -1)
      {
        Eigen::Vector3d from_parent = branches[node.id].centre - branches[par].centre;
        if (branches[node.id].dir.dot(from_parent) < 0.0)
          branches[node.id].dir *= -1.0;
      }
      for (int i = 0; i<search_size && indices(i, node.id) > -1; i++)
      {
        int child = indices(i, node.id);
        // branches[node.id].centre.setZero();
        // branches[node.id].dir = Eigen::Vector3d(0,0,1);
        // branches[child].centre = Eigen::Vector3d(1.0, 0.5, 1.0);
        // branches[child].dir = Eigen::Vector3d(1,0,0);
        
        Eigen::Vector3d to_child = branches[child].centre - branches[node.id].centre;
        // we cylinder dirs could be the opposite direction, we don't know yet
        Eigen::Vector3d dir1 = branches[node.id].dir;
        Eigen::Vector3d dir2 = branches[child].dir;
        if (dir2.dot(to_child) < 0.0)
          dir2 = -dir2;
#if 1  // Good but it needs testing!
        Eigen::Vector3d across = dir1.cross(dir2).normalized();
        Eigen::Vector3d orth1 = across.cross(dir1);
        Eigen::Vector3d orth2 = across.cross(dir2); 

        double d1 = to_child.dot(orth2)/dir1.dot(orth2);
        d1 = std::max(0.25*branches[node.id].length, d1);
        Eigen::Vector3d mid1 = branches[node.id].centre + dir1*d1;
        double d2 = to_child.dot(orth1)/dir2.dot(orth1);
        d2 = std::max(0.25*branches[child].length, d2);
        Eigen::Vector3d mid2 = branches[child].centre - dir2*d2;

        const double scale = 3.0;
        double dist = (mid1 - branches[node.id].centre).norm() + scale*(mid2 - mid1).norm() + (branches[child].centre - mid2).norm();
 //       std::cout << "d1: " << d1 << ", d2: " << d2 << ", dist: " << dist << std::endl;       
#else
        // generate square distance from 1st order approximation of the Bezier curve between branches
        double length = to_child.norm();
        double len1 = std::max(length/3.0, branches[node.id].length/4.0);
        double len2 = std::max(length/3.0, branches[child].length/4.0);
        Eigen::Vector3d mid1 = branches[node.id].centre + dir1*len1;
        Eigen::Vector3d mid2 = branches[child].centre - dir2*len2;

        const double scale = 2.0; // larger requires more straightness
        double iscale = 1.0/scale;
        Eigen::Vector3d vec1 = mid1 - branches[node.id].centre;
        Eigen::Vector3d vec2 = mid2 - mid1;
        Eigen::Vector3d vec3 = branches[child].centre - mid2;
        vec1 += dir1*(iscale - 1.0)*dir1.dot(vec1);
        vec2 += dir1*(iscale - 1.0)*dir1.dot(vec2);
        vec3 += dir1*(iscale - 1.0)*dir1.dot(vec3);
        double dist = (vec1.norm() + vec2.norm() + vec3.norm())*scale;
#endif
        double dist_sqr = dist*dist / branches[node.id].radius;

        #if defined MINIMISE_ANGLE
        const double power = 1.0;
        dist_sqr /= std::pow(std::max(0.001, to_child.dot(dir) / length), power);
        #endif
        double new_score = node.score + dist_sqr;
        if (new_score < branches[child].tree_score)
        {
					branches[child].tree_score = new_score;
					branches[child].distance_to_ground = branches[node.id].distance_to_ground + dist;
          branches[child].parent = node.id;
					closest_node.push(QueueNode(branches[child].tree_score, child));
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
    
    Eigen::Vector3d to_child = branch.centre - branches[branch.parent].centre;
    // we cylinder dirs could be the opposite direction, we don't know yet
    Eigen::Vector3d dir = branches[branch.parent].dir;
    Eigen::Vector3d dir_child = branch.dir;

    // generate square distance from 1st order approximation of the Bezier curve between branches
    double length = to_child.norm();
    double len1 = std::max(length/3.0, branches[branch.parent].length/2.0);
    double len2 = std::max(length/3.0, branch.length/2.0);
    Eigen::Vector3d mid1 = branches[branch.parent].centre + dir*len1;
    Eigen::Vector3d mid2 = branch.centre - dir_child*len2;

    Eigen::Vector3d col;
    col[0] = std::fmod(branch.tree_score,       1.0);
    col[1] = std::fmod(branch.tree_score/10.0,  1.0);
    col[2] = std::fmod(branch.tree_score/100.0, 1.0);
    
    starts.push_back(branches[branch.parent].centre);
    ends.push_back(mid1);
    colours.push_back(col);
    
    starts.push_back(mid1);
    ends.push_back(mid2);
    colours.push_back(col);
    
    starts.push_back(mid2);
    ends.push_back(branch.centre);
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
  ofs << "# tree branches file:" << std::endl;
  ofs << "x,y,z,radius" << std::endl;
  for (auto &branch: best_branches)
  {
    if (!branch.active)
      continue;
    Eigen::Vector3d base = branch.centre - branch.dir*branch.length*0.5;
    ofs << base[0] << ", " << base[1] << ", " << base[2] << ", " << branch.radius << std::endl;
  }
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
