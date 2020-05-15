// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
#include "raydebugdraw.h"
#include <nabo/nabo.h>
#include <set>

using namespace std;
using namespace Eigen;
using namespace ray;

void Cloud::save(const std::string &file_name)
{
  string name = file_name;
  if (name.substr(name.length() - 4) != ".ply")
    name += ".ply";
  writePly(name, starts, ends, times, colours);
}

bool Cloud::load(const std::string &file_name)
{
  // look first for the raycloud PLY
  if (file_name.substr(file_name.size() - 4) == ".ply")
    return loadPLY(file_name);
  if (ifstream((file_name + ".ply").c_str(), ios::in))
    return loadPLY(file_name + ".ply");

  // otherwise, look for a .laz and _traj.txt file by that name
  if (ifstream((file_name + ".laz").c_str(), ios::in) && ifstream((file_name + "_traj.txt").c_str(), ios::in))
    return loadLazTraj(file_name + ".laz", file_name + "_traj.txt");

  return false;
}

bool Cloud::load(const std::string &point_cloud, const std::string &traj_file)
{
  string name_end = point_cloud.substr(point_cloud.size() - 4);
  if (name_end == ".ply")
    readPly(point_cloud, starts, ends, times, colours);
  else if (name_end == ".laz" || name_end == ".las")
    readLas(point_cloud, ends, times, colours, 1);
  else
  {
    cout << "Error converting unknown type: " << point_cloud << endl;
    return false;
  }

  Trajectory trajectory;
  trajectory.load(traj_file);
  calculateStarts(trajectory);
  return true;
}

bool Cloud::loadPLY(const string &file)
{
  return readPly(file, starts, ends, times, colours);
}

bool Cloud::loadLazTraj(const string &laz_file, const string &traj_file)
{
  bool success = readLas(laz_file, ends, times, colours, 1);
  if (!success)
    return false;
  Trajectory trajectory;
  trajectory.load(traj_file);
  calculateStarts(trajectory);
  return true;
}

void Cloud::calculateStarts(const Trajectory &trajectory)
{
  // Aha!, problem in calculating starts when times are not ordered.
  if (trajectory.nodes.size() > 0)
  {
    int n = 1;
    starts.resize(ends.size());
    for (size_t i = 0; i < ends.size(); i++)
    {
      while ((times[i] > trajectory.nodes[n].time) && n < (int)trajectory.nodes.size() - 1) n++;
      double blend =
        (times[i] - trajectory.nodes[n - 1].time) / (trajectory.nodes[n].time - trajectory.nodes[n - 1].time);
      starts[i] =
        trajectory.nodes[n - 1].pose.position +
        (trajectory.nodes[n].pose.position - trajectory.nodes[n - 1].pose.position) * clamped(blend, 0.0, 1.0);
    }
  }
  else
    cout << "can only recalculate when a trajectory is available" << endl;
}

Vector3d Cloud::calcMinBound()
{
  Vector3d min_v(1e10, 1e10, 1e10);
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      min_v = minVector(min_v, ends[i]);
  }
  return min_v;
}

Vector3d Cloud::calcMaxBound()
{
  Vector3d max_v(-1e10, -1e10, -1e10);
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (rayBounded(i))
      max_v = maxVector(max_v, ends[i]);
  }
  return max_v;
}

void Cloud::transform(const Pose &pose, double time_delta)
{
  for (int i = 0; i < (int)starts.size(); i++)
  {
    starts[i] = pose * starts[i];
    ends[i] = pose * ends[i];
    times[i] += time_delta;
  }
}

void Cloud::removeUnboundedRays()
{
  vector<int> valids;
  for (int i = 0; i < (int)ends.size(); i++)
    if (rayBounded(i))
      valids.push_back(i);
  for (int i = 0; i < (int)valids.size(); i++)
  {
    starts[i] = starts[valids[i]];
    ends[i] = ends[valids[i]];
    times[i] = times[valids[i]];
    colours[i] = colours[valids[i]];
  }
  starts.resize(valids.size());
  ends.resize(valids.size());
  times.resize(valids.size());
  colours.resize(valids.size());
}

void Cloud::decimate(double voxel_width)
{
  vector<int> subsample = voxelSubsample(ends, voxel_width);
  for (int i = 0; i < (int)subsample.size(); i++)
  {
    int id = subsample[i];
    starts[i] = starts[id];
    ends[i] = ends[id];
    colours[i] = colours[id];
    times[i] = times[id];
  }
  starts.resize(subsample.size());
  ends.resize(subsample.size());
  colours.resize(subsample.size());
  times.resize(subsample.size());
}

void Cloud::getSurfels(int search_size, vector<Vector3d> *centroids, vector<Vector3d> *normals,
                       vector<Vector3d> *dimensions, vector<Matrix3d> *mats, MatrixXi *neighbour_indices)
{
  // simplest scheme... find 3 nearest neighbours and do cross product
  if (centroids)
    centroids->resize(ends.size());
  if (normals)
    normals->resize(ends.size());
  if (dimensions)
    dimensions->resize(ends.size());
  if (mats)
    mats->resize(ends.size());
  Nabo::NNSearchD *nns;
  vector<int> ray_ids;
  ray_ids.reserve(ends.size());
  for (unsigned int i = 0; i < ends.size(); i++)
    if (rayBounded(i))
      ray_ids.push_back(i);
  MatrixXd points_p(3, ray_ids.size());
  for (unsigned int i = 0; i < ray_ids.size(); i++) points_p.col(i) = ends[ray_ids[i]];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(search_size, ray_ids.size());
  dists2.resize(search_size, ray_ids.size());
  nns->knn(points_p, indices, dists2, search_size, 0.01, 0, 1.0);  // TODO: needs to sort here
  delete nns;

  if (neighbour_indices)
    neighbour_indices->resize(search_size, ends.size());
  for (int i = 0; i < (int)ray_ids.size(); i++)
  {
    int ii = ray_ids[i];
    if (neighbour_indices)
    {
      int j;
      for (j = 0; j < search_size && indices(j, i) > -1; j++) (*neighbour_indices)(j, ii) = ray_ids[indices(j, i)];
      if (j < search_size)
        (*neighbour_indices)(j, ii) = -1;
    }

    Vector3d centroid = ends[i];
    int num;
    for (num = 0; num < search_size && indices(num, i) > -1; num++) centroid += ends[ray_ids[indices(num, i)]];
    centroid /= (double)(num + 1);
    if (centroids)
      (*centroids)[i] = centroid;

    Matrix3d scatter = (ends[i] - centroid) * (ends[i] - centroid).transpose();
    for (int j = 0; j < num; j++)
    {
      Vector3d offset = ends[ray_ids[indices(j, i)]] - centroid;
      scatter += offset * offset.transpose();
    }
    scatter /= (double)(num + 1);

    SelfAdjointEigenSolver<Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Success);
    if (normals)
    {
      Vector3d normal = eigen_solver.eigenvectors().col(0);
      if ((ends[i] - starts[i]).dot(normal) > 0.0)
        normal = -normal;
      (*normals)[i] = normal;
    }
    if (dimensions)
    {
      Vector3d eigenvals = maxVector(Vector3d(1e-10, 1e-10, 1e-10), eigen_solver.eigenvalues());
      (*dimensions)[i] = Vector3d(sqrt(eigenvals[0]), sqrt(eigenvals[1]), sqrt(eigenvals[2]));
    }
    if (mats)
      (*mats)[i] = eigen_solver.eigenvectors();
  }
}

// starts are required to get the normal the right way around
vector<Vector3d> Cloud::generateNormals(int search_size)
{
  vector<Vector3d> normals;
  getSurfels(search_size, NULL, &normals, NULL, NULL, NULL);
  return normals;
}

// TODO: this could call getSurfels too!
void Cloud::generateEllipsoids(vector<Ellipsoid> &ellipsoids)
{
  cout << "generating " << ends.size() << " ellipsoids" << endl;
  ellipsoids.resize(ends.size());
  int search_size = 16;
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  MatrixXd points_p(3, ends.size());
  for (unsigned int i = 0; i < ends.size(); i++) points_p.col(i) = ends[i];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  MatrixXi indices;
  MatrixXd dists2;
  indices.resize(search_size, ends.size());
  dists2.resize(search_size, ends.size());
  nns->knn(points_p, indices, dists2, search_size, 0.01, 0, 1.0);
  delete nns;

  for (unsigned int i = 0; i < ends.size(); i++)
  {
    ellipsoids[i].transient = false;
    if (!rayBounded(i))
    {
      ellipsoids[i].extents.setZero();
      continue;
    }
    Matrix3d scatter;
    scatter.setZero();
    Vector3d centroid(0, 0, 0);
    double num_neighbours = 1e-10;
    for (int j = 0; j < search_size && indices(j, i) > -1; j++)
    {
      int index = indices(j, i);
      if (rayBounded(index))
      {
        centroid += ends[index];
        num_neighbours++;
      }
    }
    centroid /= num_neighbours;
    for (int j = 0; j < search_size && indices(j, i) > -1; j++)
    {
      int index = indices(j, i);
      if (rayBounded(index))
      {
        Vector3d offset = ends[index] - centroid;
        scatter += offset * offset.transpose();
      }
    }
    scatter /= num_neighbours;

    SelfAdjointEigenSolver<Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Success);

    Vector3d eigen_value = eigen_solver.eigenvalues();
    Matrix3d eigen_vector = eigen_solver.eigenvectors();

    ellipsoids[i].pos = centroid;
    double scale = 1.7;  // this scale roughly matches the dimensions of a uniformly dense ellipsoid
    eigen_value[0] = scale * sqrt(max(1e-10, eigen_value[0]));
    eigen_value[1] = scale * sqrt(max(1e-10, eigen_value[1]));
    eigen_value[2] = scale * sqrt(max(1e-10, eigen_value[2]));
    ellipsoids[i].eigen_mat.row(0) = eigen_vector.col(0) / eigen_value[0];
    ellipsoids[i].eigen_mat.row(1) = eigen_vector.col(1) / eigen_value[1];
    ellipsoids[i].eigen_mat.row(2) = eigen_vector.col(2) / eigen_value[2];
    ellipsoids[i].time = times[i];
    ellipsoids[i].setExtents(eigen_vector, eigen_value);
    ellipsoids[i].setPlanarity(eigen_value);
  }
}

void fillGrid(Grid<int> &grid, const vector<Vector3d> &starts, const vector<Vector3d> &ends)
{
  cout << "filling grid with " << ends.size() << " rays" << endl;
  // next populate the grid with these ellipsoid centres
  for (int i = 0; i < (int)ends.size(); i++)
  {
    if (!(i % 20000))
      cout << i << "/" << ends.size() << endl;
    Vector3d dir = ends[i] - starts[i];
    Vector3d dir_sign(sgn(dir[0]), sgn(dir[1]), sgn(dir[2]));
    Vector3d start = (starts[i] - grid.box_min) / grid.voxel_width;
    Vector3d end = (ends[i] - grid.box_min) / grid.voxel_width;
    Vector3i start_index(start.cast<int>());
    Vector3i end_index(end.cast<int>());
    double length_sqr = (end_index - start_index).squaredNorm();
    Vector3i index = start_index;
    for (;;)
    {
      grid.insert(index[0], index[1], index[2], i);

      if (index == end_index || (index - start_index).squaredNorm() > length_sqr)
        break;
      Vector3d mid = grid.box_min + grid.voxel_width * Vector3d(index[0] + 0.5, index[1] + 0.5, index[2] + 0.5);
      Vector3d next_boundary = mid + 0.5 * grid.voxel_width * dir_sign;
      Vector3d delta = next_boundary - starts[i];
      Vector3d d(delta[0] / dir[0], delta[1] / dir[1], delta[2] / dir[2]);
      if (d[0] < d[1] && d[0] < d[2])
        index[0] += int(dir_sign[0]);
      else if (d[1] < d[0] && d[1] < d[2])
        index[1] += int(dir_sign[1]);
      else
        index[2] += int(dir_sign[2]);
    }
  }

  grid.report();
}

void Cloud::markIntersectedEllipsoids(Grid<int> &grid, vector<bool> &transients, vector<Ellipsoid> &ellipsoids,
                                      const string &merge_type, double num_rays, bool self_transient)
{
  cout << "mark intersected ellipsoids, num_rays: " << num_rays << ", merge_type: " << merge_type << endl;
  if (DebugDraw::instance())
  {
    DebugDraw::instance()->drawCloud(ends, 1.0, 0);
  }
  int type = merge_type == "oldest" ? 0 : (merge_type == "newest" ? 1 : (merge_type == "min" ? 2 : 3));
  vector<bool> ray_tested;
  ray_tested.resize(ends.size(), false);
  int cnt = 0;
  for (auto &ellipsoid : ellipsoids)
  {
    if ((cnt++) % 20000 == 0)
      cout << cnt << "/" << ellipsoids.size() << endl;

    //    if (num_rays>0 && (ellipsoid.pos - Vector3d(2,-0.5,0.5)).norm() < 0.2)
    //      cout << "point" << endl;
    if (ellipsoid.transient)  // a previous ellipsoid could have removed the ray that represents this ellipsoid
      continue;
    if (ellipsoid.extents == Vector3d::Zero())  // unbounded rays cannot be a transient object
      continue;
    // get all the rays that overlap this ellipsoid
    Vector3d b_min = (ellipsoid.pos - ellipsoid.extents - grid.box_min) / grid.voxel_width;
    Vector3d b_max = (ellipsoid.pos + ellipsoid.extents - grid.box_min) / grid.voxel_width;
    if (b_max[0] < 0.0 || b_max[1] < 0.0 || b_max[2] < 0.0)
      continue;
    if (b_min[0] > (double)grid.dims[0] - 1.0 || b_min[1] > (double)grid.dims[1] - 1.0 ||
        b_min[2] > (double)grid.dims[2] - 1.0)
      continue;
    Vector3i bmin = maxVector(Vector3i(0, 0, 0), Vector3i(b_min.cast<int>()));
    Vector3i bmax =
      minVector(Vector3i(b_max.cast<int>()), Vector3i(grid.dims[0] - 1, grid.dims[1] - 1, grid.dims[2] - 1));

    vector<int> ray_ids;
    for (int x = bmin[0]; x <= bmax[0]; x++)
    {
      for (int y = bmin[1]; y <= bmax[1]; y++)
      {
        for (int z = bmin[2]; z <= bmax[2]; z++)
        {
          auto &list = grid.cell(x, y, z).data;
          for (auto &i : list)
          {
            if (ray_tested[i])
              continue;
            ray_tested[i] = true;
            ray_ids.push_back(i);
          }
        }
      }
    }
    for (auto &ray_id : ray_ids) ray_tested[ray_id] = false;

    double first_intersection_time = 1e10;
    double last_intersection_time = -1e10;
    int hits = 0;
    vector<int> pass_through_ids;
    for (auto &ray_id : ray_ids)
    {
      Vector3d dir = ends[ray_id] - starts[ray_id];
      // ray-ellipsoid intersection
      Vector3d to_sphere = ellipsoid.pos - starts[ray_id];
      Vector3d ray = ellipsoid.eigen_mat * dir;
      double ray_length_sqr = ray.squaredNorm();
      Vector3d to = ellipsoid.eigen_mat * to_sphere;

      double d = to.dot(ray)/ray_length_sqr;
      double dist2 = (to - ray * d).squaredNorm();

      if (dist2 > 1.0)  // misses the ellipsoid
        continue;
      double along_dist = sqrt(1.0 - dist2);
      double ray_length = sqrt(ray_length_sqr);
      d *= ray_length;
      if (ray_length < d - along_dist)  // doesn't reach the ellipsoid
        continue;

      const double pass_distance = 0.05;
      double ratio = pass_distance / dir.norm();
      bool pass_through =
        ray_length * (1.0 - ratio) > d + along_dist;  // last number requires rays to pass some way past the object
      if (pass_through)
        pass_through_ids.push_back(ray_id);
      else
      {
        hits++;
        first_intersection_time = min(first_intersection_time, times[ray_id]);
        last_intersection_time = max(last_intersection_time, times[ray_id]);
      }
    }
    size_t num_before = 0, num_after = 0;
    ellipsoid.num_rays = hits + pass_through_ids.size();
    if (num_rays == 0 || self_transient)
      ellipsoid.opacity = (double)hits / ((double)hits + (double)pass_through_ids.size());
    if (ellipsoid.num_rays == 0 || ellipsoid.opacity == 0 || num_rays == 0)
      continue;
    if (self_transient)
    {
      ellipsoid.num_gone = pass_through_ids.size();
      // now get some density stats...
      double misses = 0;
      for (auto &ray_id : pass_through_ids)
      {
        if (times[ray_id] > last_intersection_time)
          num_after++;
        else if (times[ray_id] < first_intersection_time)
          num_before++;
        else
          misses++;
      }
      double h = hits + 1e-8 - 1.0;  // subtracting 1 gives an unbiased opacity estimate
      ellipsoid.opacity = h / (h + misses);
      ellipsoid.num_gone = num_before + num_after;
    }
    else  // compare to other cloud
    {
      if (pass_through_ids.size() > 0)
      {
        if (times[pass_through_ids[0]] > ellipsoid.time)
          num_after = pass_through_ids.size();
        else
          num_before = pass_through_ids.size();
      }
    }

    double sequence_length = num_rays / ellipsoid.opacity;
    int remove_ellipsoid = false;
    if (type == 0 || type == 1)
    {
      if (double(std::max(num_before, num_after)) < sequence_length)
        continue;
      if (type == 0)                                               // oldest
        remove_ellipsoid = double(num_before) >= sequence_length;  // if false then remove numAfter rays if > seqLength
      else if (type == 1)                                          // newest
        remove_ellipsoid = double(num_after) >= sequence_length;   // if false then remove numBefore rays if > seqLength
    }
    else
    {
      // we use sum rather than max below, because it better picks out moving objects that may have some
      // pass through rays before and after the hit points.
      if (double(num_before + num_after) < sequence_length)
        continue;
      remove_ellipsoid = type == 2;  // min is remove ellipsoid, max is remove ray
    }

    if (remove_ellipsoid)
      ellipsoid.transient = true;
    else  // if we don't remove the ellipsoid then we should remove numBefore and numAfter rays if they're greater than
          // sequence length
    {
      double d = 0.0;
      for (int j = 0; j < (int)pass_through_ids.size(); j++)
      {
        d += ellipsoid.opacity;
        if (d >= 1.0)
          d--;
        else
          continue;
        int i = pass_through_ids[j];
        if (!self_transient || times[i] < first_intersection_time || times[i] > last_intersection_time)
        {
          // remove ray i
          transients[i] = true;
        }
      }
    }
  }
}

static double voxel_width = 0.25;

void Cloud::findTransients(Cloud &transient, Cloud &fixed, const string &merge_type, double num_rays, bool colour_cloud)
{
  cout << "find transients" << endl;
  Vector3d box_min(1e10, 1e10, 1e10);
  Vector3d box_max(-1e10, -1e10, -1e10);
  for (auto &end : ends)
  {
    box_min = minVector(box_min, end);
    box_max = maxVector(box_max, end);
  }

  vector<Ellipsoid> ellipsoids;
  generateEllipsoids(ellipsoids);

  Grid<int> grid(box_min, box_max, voxel_width);
  fillGrid(grid, starts, ends);

  vector<bool> transients;
  transients.resize(ends.size(), false);
  // now walk every ray through the grid and mark if transient
  markIntersectedEllipsoids(grid, transients, ellipsoids, merge_type, num_rays, true);

  // Lastly, generate the new ray clouds from this sphere information
  for (int i = 0; i < (int)ellipsoids.size(); i++)
  {
    RGBA col = colours[i];
    if (colour_cloud)
    {
      col.red = (uint8_t)((1.0 - ellipsoids[i].planarity) * 255.0);
      col.blue = (uint8_t)(ellipsoids[i].opacity * 255.0);
      col.green = (uint8_t)((double)ellipsoids[i].num_gone / ((double)ellipsoids[i].num_gone + 10.0) * 255.0);
    }
    if (ellipsoids[i].transient || transients[i])
    {
      transient.starts.push_back(starts[i]);
      transient.ends.push_back(ends[i]);
      transient.times.push_back(times[i]);
      transient.colours.push_back(col);
    }
    else
    {
      fixed.starts.push_back(starts[i]);
      fixed.ends.push_back(ends[i]);
      fixed.times.push_back(times[i]);
      fixed.colours.push_back(col);
    }
  }
}

void Cloud::combine(std::vector<Cloud> &clouds, Cloud &differences, const string &merge_type, double num_rays)
{
  vector<Grid<int>> grids(clouds.size());
  for (int c = 0; c < (int)clouds.size(); c++)
  {
    Vector3d box_min(1e10, 1e10, 1e10);
    Vector3d box_max(-1e10, -1e10, -1e10);
    for (auto &end : clouds[c].ends)
    {
      box_min = minVector(box_min, end);
      box_max = maxVector(box_max, end);
    }
    grids[c].init(box_min, box_max, voxel_width);
    fillGrid(grids[c], clouds[c].starts, clouds[c].ends);
  }

  vector<vector<bool>> transients(clouds.size());
  for (int c = 0; c < (int)clouds.size(); c++) transients[c].resize(clouds[c].ends.size(), false);
  // now for each cloud, look for other clouds that penetrate it
  for (int c = 0; c < (int)clouds.size(); c++)
  {
    vector<Ellipsoid> ellipsoids;
    clouds[c].generateEllipsoids(ellipsoids);
    // just set opacity
    clouds[c].markIntersectedEllipsoids(grids[c], transients[c], ellipsoids, merge_type, 0, false);

    for (int d = 0; d < (int)clouds.size(); d++)
    {
      if (d == c)
        continue;
      // use ellipsoid opacity to set transient flag true on transients
      clouds[d].markIntersectedEllipsoids(grids[d], transients[d], ellipsoids, merge_type, num_rays, false);
    }

    for (int i = 0; i < (int)clouds[c].ends.size(); i++)
      if (ellipsoids[i].transient)
        transients[c][i] = true;  // HACK, we need a better way to signal this!
  }
  for (int c = 0; c < (int)clouds.size(); c++)
  {
    auto &cloud = clouds[c];
    int t = 0;
    int f = 0;
    for (int i = 0; i < (int)cloud.ends.size(); i++)
    {
      if (transients[c][i])
      {
        t++;
        differences.starts.push_back(cloud.starts[i]);
        differences.ends.push_back(cloud.ends[i]);
        differences.times.push_back(cloud.times[i]);
        differences.colours.push_back(cloud.colours[i]);
      }
      else
      {
        f++;
        starts.push_back(cloud.starts[i]);
        ends.push_back(cloud.ends[i]);
        times.push_back(cloud.times[i]);
        colours.push_back(cloud.colours[i]);
      }
    }
    cout << t << " transients, " << f << " fixed rays." << endl;
  }
}

void Cloud::split(Cloud &cloud1, Cloud &cloud2, function<bool(int i)> fptr)
{
  for (int i = 0; i < (int)ends.size(); i++)
  {
    Cloud &cloud = fptr(i) ? cloud2 : cloud1;
    cloud.starts.push_back(starts[i]);
    cloud.ends.push_back(ends[i]);
    cloud.times.push_back(times[i]);
    cloud.colours.push_back(colours[i]);
  }
}
