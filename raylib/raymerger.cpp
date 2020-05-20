#include "raymerger.h"
#include <limits>
#include "raytrajectory.h"
#include "rayply.h"
#include "raylaz.h"
#include "raydebugdraw.h"
#include <nabo/nabo.h>

namespace
{
typedef Eigen::Matrix<double, 6, 1> Vector6i;
struct Vector6iLess
{
  inline bool operator()(const Vector6i &a, const Vector6i &b) const
  {
    if (a[0] != b[0])
      return a[0] < b[0];
    if (a[1] != b[1])
      return a[1] < b[1];
    if (a[2] != b[2])
      return a[2] < b[2];
    if (a[3] != b[3])
      return a[3] < b[3];
    if (a[4] != b[4])
      return a[4] < b[4];
    return a[5] < b[5];
  }
};
const double test_width = 0.01; // allows a minor variation when checking for similarity of rays

void rayLookup(const ray::Cloud *cloud, std::set<Vector6i, Vector6iLess> &ray_lookup)
{
  for (size_t i = 0; i<cloud->ends.size(); i++)
  {
    const Eigen::Vector3d &point = cloud->ends[i];
    const Eigen::Vector3d &start = cloud->starts[i];
    Vector6i ray;
    for (int j = 0; j<3; j++)
    {
      ray[j]   = int(floor(start[j] / test_width)); 
      ray[3+j] = int(floor(point[j] / test_width));
    }
    if (ray_lookup.find(ray) == ray_lookup.end())
      ray_lookup.insert(ray);
  }
}

struct Ellipsoid
{
  Eigen::Vector3d pos;
  Eigen::Matrix3d eigen_mat; // each row is a scaled eigenvector
  double time;
  Eigen::Vector3d extents;
  double opacity;
  double planarity;
  size_t num_rays;
  size_t num_gone;
  inline void setExtents(const Eigen::Matrix3d &vecs, const Eigen::Vector3d &vals)
  {
    // This is approximate (slightly larger than minimal bounds), but
    // an exact bounding box is most likely non-analytic, and expensive to compute
    double max_rr = std::max(vals[0], std::max(vals[1], vals[2]));
    const Eigen::Vector3d &x = vecs.col(0);
    const Eigen::Vector3d &y = vecs.col(1);
    const Eigen::Vector3d &z = vecs.col(2);
    extents[0] = std::min(max_rr, abs(x[0]) * vals[0] + abs(y[0]) * vals[1] + abs(z[0]) * vals[2]);
    extents[1] = std::min(max_rr, abs(x[1]) * vals[0] + abs(y[1]) * vals[1] + abs(z[1]) * vals[2]);
    extents[2] = std::min(max_rr, abs(x[2]) * vals[0] + abs(y[2]) * vals[1] + abs(z[2]) * vals[2]);
  }
  void setPlanarity(const Eigen::Vector3d &vals) { planarity = (vals[1] - vals[0]) / vals[1]; }
  bool transient;
};

// Convert the cloud into a list of ellipsoids, which represent a volume around each cloud point, 
// shaped by the distribution of its neighbouring points.
void generateEllipsoids(const ray::Cloud &cloud, std::vector<Ellipsoid> &ellipsoids)
{
  std::cout << "generating " << cloud.ends.size() << " ellipsoids" << std::endl;
  ellipsoids.resize(cloud.ends.size());
  int search_size = 16;
  Nabo::NNSearchD *nns;
  Nabo::Parameters params("bucketSize", 8);
  Eigen::MatrixXd points_p(3, cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++) 
    points_p.col(i) = cloud.ends[i];
  nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);

  // Run the search
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, cloud.ends.size());
  dists2.resize(search_size, cloud.ends.size());
  nns->knn(points_p, indices, dists2, search_size, 0.01, 0, 1.0);
  delete nns;

  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    ellipsoids[i].transient = false;
    ellipsoids[i].opacity = 1.0;
    if (!cloud.rayBounded(i))
    {
      ellipsoids[i].extents.setZero();
      continue;
    }
    Eigen::Matrix3d scatter;
    scatter.setZero();
    Eigen::Vector3d centroid(0, 0, 0);
    double num_neighbours = 0;
    for (int j = 0; j < search_size && indices(j, i) > -1; j++)
    {
      int index = indices(j, i);
      if (cloud.rayBounded(index))
      {
        centroid += cloud.ends[index];
        num_neighbours++;
      }
    }
    if (num_neighbours < 4)
    {
      ellipsoids[i].extents.setZero();
      continue;
    }
    centroid /= num_neighbours;
    for (int j = 0; j < search_size && indices(j, i) > -1; j++)
    {
      int index = indices(j, i);
      if (cloud.rayBounded(index))
      {
        Eigen::Vector3d offset = cloud.ends[index] - centroid;
        scatter += offset * offset.transpose();
      }
    }
    scatter /= num_neighbours;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());
    ASSERT(eigen_solver.info() == Success);

    Eigen::Vector3d eigen_value = eigen_solver.eigenvalues();
    Eigen::Matrix3d eigen_vector = eigen_solver.eigenvectors();

    ellipsoids[i].pos = centroid;
    double scale = 1.7;  // this scale roughly matches the dimensions of a uniformly dense ellipsoid
    eigen_value[0] = scale * std::sqrt(std::max(1e-10, eigen_value[0]));
    eigen_value[1] = scale * std::sqrt(std::max(1e-10, eigen_value[1]));
    eigen_value[2] = scale * std::sqrt(std::max(1e-10, eigen_value[2]));
    ellipsoids[i].eigen_mat.row(0) = eigen_vector.col(0) / eigen_value[0];
    ellipsoids[i].eigen_mat.row(1) = eigen_vector.col(1) / eigen_value[1];
    ellipsoids[i].eigen_mat.row(2) = eigen_vector.col(2) / eigen_value[2];
    ellipsoids[i].time = cloud.times[i];
    ellipsoids[i].setExtents(eigen_vector, eigen_value);
    ellipsoids[i].setPlanarity(eigen_value);
  }
}

// Flag each ellipsoid that has been intersected by a ray in the grid structure, for removal. 
// Depending on the merge_type, we may instead flag the ray for removal instead (using the transients vector)
void markIntersectedEllipsoids(const ray::Cloud &cloud, ray::Grid<int> &grid, std::vector<bool> &transients, std::vector<Ellipsoid> &ellipsoids,
                               const std::string &merge_type, double num_rays, bool self_transient)
{
  std::cout << "mark intersected ellipsoids, num_rays: " << num_rays << ", merge_type: " << merge_type << std::endl;
  if (ray::DebugDraw::instance())
  {
    ray::DebugDraw::instance()->drawCloud(cloud.ends, 1.0, 0);
  }
  int type = merge_type == "oldest" ? 0 : (merge_type == "newest" ? 1 : (merge_type == "min" ? 2 : 3));
  std::vector<bool> ray_tested;
  ray_tested.resize(cloud.ends.size(), false);
  int cnt = 0;
  for (auto &ellipsoid : ellipsoids)
  {
    if ((cnt++) % 20000 == 0)
      std::cout << cnt << "/" << ellipsoids.size() << std::endl;

    if (ellipsoid.transient)  // a previous ellipsoid could have removed the ray that represents this ellipsoid
      continue;
    if (ellipsoid.extents == Eigen::Vector3d::Zero())  // unbounded rays cannot be a transient object
      continue;
    // get all the rays that overlap this ellipsoid
    Eigen::Vector3d b_min = (ellipsoid.pos - ellipsoid.extents - grid.box_min) / grid.voxel_width;
    Eigen::Vector3d b_max = (ellipsoid.pos + ellipsoid.extents - grid.box_min) / grid.voxel_width;
    if (b_max[0] < 0.0 || b_max[1] < 0.0 || b_max[2] < 0.0)
      continue;
    if (b_min[0] >= (double)grid.dims[0] || b_min[1] >= (double)grid.dims[1] ||
        b_min[2] >= (double)grid.dims[2])
      continue;
    Eigen::Vector3i bmin = ray::maxVector(Eigen::Vector3i(0, 0, 0), Eigen::Vector3i(b_min.cast<int>()));
    Eigen::Vector3i bmax =
      ray::minVector(Eigen::Vector3i(b_max.cast<int>()), Eigen::Vector3i(grid.dims[0] - 1, grid.dims[1] - 1, grid.dims[2] - 1));

    std::vector<int> ray_ids;
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
    for (auto &ray_id : ray_ids) 
      ray_tested[ray_id] = false;

    double first_intersection_time = 1e10;
    double last_intersection_time = -1e10;
    int hits = 0;
    std::vector<int> pass_through_ids;
    for (auto &ray_id : ray_ids)
    {
      Eigen::Vector3d dir = cloud.ends[ray_id] - cloud.starts[ray_id];
      // ray-ellipsoid intersection
      Eigen::Vector3d to_sphere = ellipsoid.pos - cloud.starts[ray_id];
      Eigen::Vector3d ray = ellipsoid.eigen_mat * dir;
      double ray_length_sqr = ray.squaredNorm();
      Eigen::Vector3d to = ellipsoid.eigen_mat * to_sphere;

      double d = to.dot(ray)/ray_length_sqr;
      double dist2 = (to - ray * d).squaredNorm();

      if (dist2 > 1.0)  // misses the ellipsoid
        continue;
      double along_dist = std::sqrt(1.0 - dist2);
      double ray_length = std::sqrt(ray_length_sqr);
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
        first_intersection_time = std::min(first_intersection_time, cloud.times[ray_id]);
        last_intersection_time = std::max(last_intersection_time, cloud.times[ray_id]);
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
        if (cloud.times[ray_id] > last_intersection_time)
          num_after++;
        else if (cloud.times[ray_id] < first_intersection_time)
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
        if (cloud.times[pass_through_ids[0]] > ellipsoid.time)
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
      if (double(num_before + num_after) < sequence_length)  // TODO: even a tiny bit of translucency will make a single ray not enough
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
        if (!self_transient || cloud.times[i] < first_intersection_time || cloud.times[i] > last_intersection_time)
        {
          // remove ray i
          transients[i] = true;
        }
      }
    }
  }
}

// Fill the grid acceleration structure with the rays in a ray cloud
void fillGrid(ray::Grid<int> &grid, const ray::Cloud &cloud)
{
  std::cout << "filling grid with " << cloud.ends.size() << " rays" << std::endl;
  // next populate the grid with these ellipsoid centres
  for (int i = 0; i < (int)cloud.ends.size(); i++)
  {
    if (!(i % 20000))
      std::cout << i << "/" << cloud.ends.size() << std::endl;
    Eigen::Vector3d dir = cloud.ends[i] - cloud.starts[i];
    Eigen::Vector3d dir_sign(ray::sgn(dir[0]), ray::sgn(dir[1]), ray::sgn(dir[2]));
    Eigen::Vector3d start = (cloud.starts[i] - grid.box_min) / grid.voxel_width;
    Eigen::Vector3d end = (cloud.ends[i] - grid.box_min) / grid.voxel_width;
    Eigen::Vector3i start_index((int)floor(start[0]), (int)floor(start[1]), (int)floor(start[2]));
    Eigen::Vector3i end_index((int)floor(end[0]), (int)floor(end[1]), (int)floor(end[2]));
    double length_sqr = (end_index - start_index).squaredNorm();
    Eigen::Vector3i index = start_index;
    for (;;)
    {
      grid.insert(index[0], index[1], index[2], i);

      if (index == end_index || (index - start_index).squaredNorm() > length_sqr)
        break;
      Eigen::Vector3d mid = grid.box_min + grid.voxel_width * Eigen::Vector3d(index[0] + 0.5, index[1] + 0.5, index[2] + 0.5);
      Eigen::Vector3d next_boundary = mid + 0.5 * grid.voxel_width * dir_sign;
      Eigen::Vector3d delta = next_boundary - cloud.starts[i];
      Eigen::Vector3d d(delta[0] / dir[0], delta[1] / dir[1], delta[2] / dir[2]);
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

}

void ray::Merger::mergeSingleCloud(const ray::Cloud &cloud, const std::string &merge_type, double num_rays, bool colour_cloud)
{
  const double voxel_width = 4.0 * cloud.estimatePointSpacing();
  std::cout << "find transients" << std::endl;

  std::vector<Ellipsoid> ellipsoids;
  generateEllipsoids(cloud, ellipsoids);

  ray::Grid<int> grid(cloud.calcMinBound(), cloud.calcMaxBound(), voxel_width);
  fillGrid(grid, cloud);

  std::vector<bool> transients;
  transients.resize(cloud.ends.size(), false);
  // now walk every ray through the grid and mark if transient
  markIntersectedEllipsoids(cloud, grid, transients, ellipsoids, merge_type, num_rays, true);

  // Lastly, generate the new ray clouds from this sphere information
  for (int i = 0; i < (int)ellipsoids.size(); i++)
  {
    ray::RGBA col = cloud.colours[i];
    if (colour_cloud)
    {
      col.red = (uint8_t)((1.0 - ellipsoids[i].planarity) * 255.0);
      col.blue = (uint8_t)(ellipsoids[i].opacity * 255.0);
      col.green = (uint8_t)((double)ellipsoids[i].num_gone / ((double)ellipsoids[i].num_gone + 10.0) * 255.0);
    }
    if (ellipsoids[i].transient || transients[i])
    {
      differences_.starts.push_back(cloud.starts[i]);
      differences_.ends.push_back(cloud.ends[i]);
      differences_.times.push_back(cloud.times[i]);
      differences_.colours.push_back(col);
    }
    else
    {
      result_.starts.push_back(cloud.starts[i]);
      result_.ends.push_back(cloud.ends[i]);
      result_.times.push_back(cloud.times[i]);
      result_.colours.push_back(col);
    }
  }
}

void ray::Merger::mergeMultipleClouds(std::vector<ray::Cloud> &clouds, const std::string &merge_type, double num_rays)
{
  std::vector<ray::Grid<int>> grids(clouds.size());
  for (int c = 0; c < (int)clouds.size(); c++)
  {
    grids[c].init(clouds[c].calcMinBound(), clouds[c].calcMaxBound(), 4.0*clouds[c].estimatePointSpacing());
    fillGrid(grids[c], clouds[c]);
  }

  std::vector<std::vector<bool>> transients(clouds.size());
  for (int c = 0; c < (int)clouds.size(); c++) transients[c].resize(clouds[c].ends.size(), false);
  // now for each cloud, look for other clouds that penetrate it
  for (int c = 0; c < (int)clouds.size(); c++)
  {
    std::vector<Ellipsoid> ellipsoids;
    generateEllipsoids(clouds[c], ellipsoids);
    // just set opacity
    markIntersectedEllipsoids(clouds[c], grids[c], transients[c], ellipsoids, merge_type, 0, false);

    for (int d = 0; d < (int)clouds.size(); d++)
    {
      if (d == c)
        continue;
      // use ellipsoid opacity to set transient flag true on transients
      markIntersectedEllipsoids(clouds[d], grids[d], transients[d], ellipsoids, merge_type, num_rays, false);
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
        differences_.starts.push_back(cloud.starts[i]);
        differences_.ends.push_back(cloud.ends[i]);
        differences_.times.push_back(cloud.times[i]);
        differences_.colours.push_back(cloud.colours[i]);
      }
      else
      {
        f++;
        result_.starts.push_back(cloud.starts[i]);
        result_.ends.push_back(cloud.ends[i]);
        result_.times.push_back(cloud.times[i]);
        result_.colours.push_back(cloud.colours[i]);
      }
    }
    std::cout << t << " transients, " << f << " fixed rays." << std::endl;
  }
}

void ray::Merger::mergeThreeWay(const ray::Cloud &base_cloud, ray::Cloud &cloud1, ray::Cloud &cloud2, const std::string &merge_type, double num_rays)
{
  // The 3-way merge is similar to those performed on text files for version control systems. It attempts to apply the 
  // changes in both cloud 1 and cloud2 (compared to base_cloud). When there is a conflict (different changes in the same 
  // location) it resolves that according to the selected merge_type.
  // unlike with text, a change requires a small threshold, since positions are floating point values. In our case, we 
  // define a ray as unchanged when the start and end points are within the same small voxel as they were in base_cloud.
  // so the threshold is test_width.

  // generate quick lookup for the existance of a particular (quantised) ray
  ray::Cloud *clouds[2] = {&cloud1, &cloud2};
  std::set<Vector6i, Vector6iLess> base_ray_lookup;
  rayLookup(&base_cloud, base_ray_lookup);
  std::set<Vector6i, Vector6iLess> ray_lookups[2];
  for (int c = 0; c<2; c++)
    rayLookup(clouds[c], ray_lookups[c]);

  std::cout << "set size " << ray_lookups[0].size() << ", " << ray_lookups[1].size() << ", " << base_ray_lookup.size() << std::endl;

  // now remove all similar rays to base_cloud and put them in the final cloud:
  int preferred_cloud = clouds[0]->times[0] > clouds[1]->times[0] ? 0 : 1;
  int u = 0;
  for (int c = 0; c < 2; c++)
  {
    ray::Cloud &cloud = *clouds[c];
    for (int i = 0; i<(int)cloud.ends.size(); i++)
    {
      Eigen::Vector3d &point = cloud.ends[i];
      Eigen::Vector3d &start = cloud.starts[i];
      Vector6i ray;
      for (int j = 0; j<3; j++)
      {
        ray[j]   = int(std::floor(start[j] / test_width)); 
        ray[3+j] = int(std::floor(point[j] / test_width));
      }
      int other = 1-c;
      // if the ray is in cloud1 and cloud2 there is no contention, so add the ray to the result
      if (ray_lookups[other].find(ray) != ray_lookups[other].end()) 
      {
        if (c == preferred_cloud)
        {
          result_.starts.push_back(start);
          result_.ends.push_back(point);
          result_.times.push_back(cloud.times[i]);
          result_.colours.push_back(cloud.colours[i]);
          u++;
        }
      }
      // we want to run the combine (which revolves conflicts) on only the changed parts
      // so we want to keep only the changes for cloud[0] and cloud[1]...
      // which means removing rays that aren't changed: 
      if (base_ray_lookup.find(ray) != base_ray_lookup.end()) 
      {
        cloud.starts[i] = cloud.starts.back(); cloud.starts.pop_back();
        cloud.ends[i] = cloud.ends.back(); cloud.ends.pop_back();
        cloud.times[i] = cloud.times.back(); cloud.times.pop_back();
        cloud.colours[i] = cloud.colours.back(); cloud.colours.pop_back();
        i--;
      }
    }
  }
  std::cout << u << " unaltered rays have been moved into combined cloud" << std::endl;
  std::cout << clouds[0]->ends.size() << " and " << clouds[1]->ends.size() << " rays to combine, that are different" << std::endl;
#if defined VERBOSE_MERGE
  this->save("common_rays.ply");
  clouds[0]->save("changes_0.ply");
  clouds[1]->save("changes_1.ply");
#endif
  // 'all' means keep all the changes, so we simply concatenate these rays into the result
  if (merge_type == "all")
  {
    for (int c = 0; c<2; c++)
    {
      result_.starts.insert(result_.starts.end(), clouds[c]->starts.begin(), clouds[c]->starts.end());
      result_.ends.insert(result_.ends.end(), clouds[c]->ends.begin(), clouds[c]->ends.end());
      result_.times.insert(result_.times.end(), clouds[c]->times.begin(), clouds[c]->times.end());
      result_.colours.insert(result_.colours.end(), clouds[c]->colours.begin(), clouds[c]->colours.end());
    }
    return;
  }
  // otherwise we run combine on the altered clouds
  // first, grid the rays for fast lookup
  ray::Grid<int> grids[2];
  for (int c = 0; c < 2; c++)
  {
    grids[c].init(clouds[c]->calcMinBound(), clouds[c]->calcMaxBound(), 4.0*clouds[c]->estimatePointSpacing());
    fillGrid(grids[c], *clouds[c]);
  }  

  std::vector<bool> transients[2];
  for (int c = 0; c < 2; c++) 
    transients[c].resize(clouds[c]->ends.size(), false);
  // now for each cloud, represent the end points as ellipsoids, and ray cast the other cloud's rays against it
  for (int c = 0; c < 2; c++)
  {
    if (clouds[c]->ends.size()==0)
      continue;
    std::vector<Ellipsoid> ellipsoids;
    generateEllipsoids(*clouds[c], ellipsoids);

    // just set opacity
    markIntersectedEllipsoids(*clouds[c], grids[c], transients[c], ellipsoids, merge_type, 0, false);

    int d = 1 - c;
    // use ellipsoid opacity to set transient flag true on transients (intersected ellipsoids)
    markIntersectedEllipsoids(*clouds[d], grids[d], transients[d], ellipsoids, merge_type, num_rays, false);

    for (int i = 0; i < (int)ellipsoids.size(); i++)
      if (ellipsoids[i].transient)
        transients[c][i] = true;  
  }
  for (int c = 0; c < 2; c++)
  {
    auto &cloud = *clouds[c];
    int t = 0;
    int f = 0;
    for (int i = 0; i < (int)transients[c].size(); i++)
    {
      if (!transients[c][i])
      {
        f++;
        result_.starts.push_back(cloud.starts[i]);
        result_.ends.push_back(cloud.ends[i]);
        result_.times.push_back(cloud.times[i]);
        result_.colours.push_back(cloud.colours[i]);
      }
      else
        t++; // we aren't storing the differences. No current demand for this.
    }
    std::cout << t << " transients, " << f << " fixed rays." << std::endl;
  }
}
