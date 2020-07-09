// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raybuildinggen.h"

namespace ray
{
class Cuboid
{
public:
  Cuboid(const Eigen::Vector3d &min_b, const Eigen::Vector3d &max_b)
  {
    min_bound = min_b;
    max_bound = max_b;
  }
  Cuboid(){}
  Eigen::Vector3d min_bound, max_bound;

  bool rayIntersectBox(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth) const
  {
    double max_near_d = 0;
    double min_far_d = std::numeric_limits<double>::max();
    Eigen::Vector3d centre = (min_bound + max_bound) / 2.0;
    Eigen::Vector3d extent = (max_bound - min_bound) / 2.0;
    Eigen::Vector3d to_centre = centre - start;
    for (int ax = 0; ax < 3; ax++)
    {
      double s = dir[ax] > 0.0 ? 1.0 : -1.0;
      double near_d = (to_centre[ax] - s * extent[ax]) / dir[ax];
      double far_d = (to_centre[ax] + s * extent[ax]) / dir[ax];

      max_near_d = std::max(max_near_d, near_d);
      min_far_d = std::min(min_far_d, far_d);
    }
    if (max_near_d > 0.0 && max_near_d < depth && max_near_d < min_far_d)
    {
      depth = max_near_d;
      return true;
    }
    return false;
  }

  bool rayIntersectNegativeBox(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth) const
  {
    double max_near_d = 0;
    double min_far_d = 1e10;
    Eigen::Vector3d centre = (min_bound + max_bound) / 2.0;
    Eigen::Vector3d extent = (max_bound - min_bound) / 2.0;
    Eigen::Vector3d to_centre = centre - start;
    for (int ax = 0; ax < 3; ax++)
    {
      double s = dir[ax] > 0.0 ? 1.0 : -1.0;
      double near_d = (to_centre[ax] - s * extent[ax]) / dir[ax];
      double far_d = (to_centre[ax] + s * extent[ax]) / dir[ax];

      max_near_d = std::max(max_near_d, near_d);
      min_far_d = std::min(min_far_d, far_d);
    }
    if (max_near_d < min_far_d && min_far_d < depth)
    {
      depth = min_far_d;
      return true;
    }
    return false;
  }

  bool intersects(const Eigen::Vector3d &pos) const
  {
    return pos[0] > min_bound[0] && pos[1] > min_bound[1] && pos[2] > min_bound[2] && pos[0] < max_bound[0] &&
           pos[1] < max_bound[1] && pos[2] < max_bound[2];
  }

  bool overlaps(const Cuboid &other) const
  {
    bool outside = other.min_bound[0] > max_bound[0] || other.min_bound[1] > max_bound[1] || other.min_bound[2] > max_bound[2] ||
                   other.max_bound[0] < min_bound[0] || other.max_bound[1] < min_bound[1] || other.max_bound[2] < min_bound[2];
    return !outside;
  }
};

static double table_density = 0; // items per square metre
static double cupboard_density = 0; // items per metre along wall

void split(const Cuboid &cuboid, std::vector<Cuboid> &cuboids, std::vector<Cuboid> &furniture)
{
  const double room_scales[3] = {5.0, 5.0, 2.7};
  const double door_height = 2.0;
  const double door_width = 0.5;
  const double wall_width = 0.2;
  const double floor_width = 0.5;
  const double distinct_floor_likelihood = 0.55; // high has separated floors, low is a mixture of floor heights
  Eigen::Vector3d ext = cuboid.max_bound - cuboid.min_bound; 
  double sum = ext[0] + ext[1] + ext[2];
  double r = random(0.0, sum);
  int ax = r < ext[0] ? 0 : r<ext[0] + ext[1] ? 1 : 2;
  if (ext[2] > 2.0*room_scales[2] && random(0.0,1.0) < distinct_floor_likelihood)
    ax = 2.0;
 // int ax = ext[2] > ext[1] ? (ext[2] > ext[0] ? 2 : 0) : (ext[1] > ext[0] ? 1 : 0);
  double length = ext[ax]; 
 // std::cout << "room size " << (cuboid.max_bound - cuboid.min_bound).transpose() << " splitting on " << ax << " axis." << std::endl; 
  if (length > 2.0*room_scales[ax])
  {
    double split_point = cuboid.min_bound[ax] + random(room_scales[ax], length - room_scales[ax]);
    Cuboid rooms[2] = {cuboid, cuboid};
    rooms[0].max_bound[ax] = split_point;
    rooms[1].min_bound[ax] = split_point;

    if (ax == 2) // stair hole
    {
      // just a hole for now
      Cuboid hole;
      Eigen::Vector3d hole_pos = cuboid.min_bound + ext * random(0.25, 0.75);
      hole_pos[ax] = split_point;
      hole.min_bound = hole_pos - Eigen::Vector3d(1.5*door_width, 2.5*door_width, door_width);
      hole.max_bound = hole_pos + Eigen::Vector3d(1.5*door_width, 2.5*door_width, door_width);
      cuboids.push_back(hole);
    }
    else
    { 
      Cuboid door;
      Eigen::Vector3d door_pos = cuboid.min_bound + ext * random(0.25, 0.75);
      door_pos[ax] = split_point;
      door.min_bound = door_pos - Eigen::Vector3d(door_width, door_width, 0.0);
      door.max_bound = door_pos + Eigen::Vector3d(door_width, door_width, 0.0);
      door.min_bound[2] = cuboid.min_bound[2];
      door.max_bound[2] = cuboid.min_bound[2] + door_height;
      cuboids.push_back(door);
    }
    split(rooms[0], cuboids, furniture);
    split(rooms[1], cuboids, furniture);
  }
  else
  {
    Cuboid room = cuboid;
    room.min_bound += 0.5*Eigen::Vector3d(wall_width, wall_width, 0);
    room.max_bound -= 0.5*Eigen::Vector3d(wall_width, wall_width, 2.0*floor_width);

    // now lets add some furniture...
    // tables:
    double floor_area = (room.max_bound[0]-room.min_bound[0]) * (room.max_bound[1]-room.min_bound[1]);
    int num_tables = (int)(floor_area * table_density);
    for (int i = 0; i<num_tables; i++)
    {
      double table_width = 1.0;
      double table_length = 1.5;
      const double table_height = 1.1;
      if (std::rand()%2)
        std::swap(table_width, table_length);
      Eigen::Vector3d table_rad = 0.5 * Eigen::Vector3d(table_width, table_length, table_height);
      Eigen::Vector3d table_pos(random(room.min_bound[0]+table_width/2.0, room.max_bound[0]-table_width/2.0),
                                random(room.min_bound[1]+table_length/2.0, room.max_bound[1]-table_length/2.0),
                                room.min_bound[2] + table_height/2.0);
      Cuboid table;
      table.min_bound = table_pos - table_rad;
      table.max_bound = table_pos + table_rad;
      bool overlaps = false;
      for (auto &other: furniture)
        if (other.overlaps(table))
          overlaps = true;
      for (auto &other: cuboids)
        if (other.overlaps(table))
          overlaps = true;
      if (!overlaps)
        furniture.push_back(table);
    }
    // cupboards:
    double wall_length = 2.0*((room.max_bound[0]-room.min_bound[0]) + (room.max_bound[1]-room.min_bound[1]));
    int num_cupboards = (int)(wall_length * cupboard_density);
    for (int i = 0; i<num_cupboards; i++)
    {
      double cupboard_width = 0.6;
      double cupboard_length = 2.0;
      double cupboard_height = random(1.0, 3.0);
      int wall_id = std::rand()%2;
      int side = std::rand()%2;
      if (wall_id == 1)
        std::swap(cupboard_width, cupboard_length);
      Eigen::Vector3d cupboard_rad = 0.5 * Eigen::Vector3d(cupboard_width, cupboard_length, cupboard_height);
      Eigen::Vector3d cupboard_pos;
      cupboard_pos[2] = room.min_bound[2] + 0.5*cupboard_height;
      if (wall_id == 0)
      {
        cupboard_pos[0] = side ? room.min_bound[0] + 0.5*cupboard_width : room.max_bound[0] - 0.5*cupboard_width;
        cupboard_pos[1] = random(room.min_bound[1]+0.5*cupboard_length, room.max_bound[1]-0.5*cupboard_length);
      }
      else
      {
        cupboard_pos[0] = random(room.min_bound[0]+0.5*cupboard_width, room.max_bound[0]-0.5*cupboard_width);
        cupboard_pos[1] = side ? room.min_bound[1] + 0.5*cupboard_length : room.max_bound[1] - 0.5*cupboard_length;
      }
      Cuboid cupboard;
      cupboard.min_bound = cupboard_pos - cupboard_rad;
      cupboard.max_bound = cupboard_pos + cupboard_rad;
      bool overlaps = false;
      for (auto &other: furniture)
        if (other.overlaps(cupboard))
          overlaps = true;
      for (auto &other: cuboids)
        if (other.overlaps(cupboard))
          overlaps = true;
      if (!overlaps)
        furniture.push_back(cupboard);
    }
    cuboids.push_back(room);
  }
}

void buildPath(const Eigen::Vector3d &start_pos, int id, int last_id, const std::vector<Cuboid> &cuboids,
               std::vector<Eigen::Vector3d> &points, double density, std::vector<int> &visited)
{
  visited[id] = 1;
  Eigen::Vector3d last_pos = start_pos;
 // std::cout << "num cuboids: " << cuboids.size() << std::endl;
  for (int i = 0; i<(int)cuboids.size(); i++)
  {
    if (i == id || i==last_id || visited[i])
      continue;
    if (cuboids[i].overlaps(cuboids[id]))
    {
      Eigen::Vector3d pos = (cuboids[i].max_bound + cuboids[i].min_bound)/2.0;

      int num = (int)((pos - last_pos).norm() * density);
      for (int j = 0; j<num; j++)
        points.push_back(last_pos + (pos - last_pos)*(double)j/(double)num);
      buildPath(pos, i, id, cuboids, points, density, visited);
      last_pos = pos;
    }
  }
  if (last_pos == start_pos) // no outlets found
  {
    Eigen::Vector3d pos = (cuboids[id].max_bound + cuboids[id].min_bound)/2.0;
    int num = (int)((pos - last_pos).norm() * density);
    for (int j = 0; j<num; j++)
      points.push_back(last_pos + (pos - last_pos)*(double)j/(double)num);
    last_pos = pos;
  }
  int num = (int)((start_pos - last_pos).norm() * density);
  for (int j = 0; j<num; j++)
    points.push_back(last_pos + (start_pos - last_pos)*(double)j/(double)num);
}

// A room with a door, window, table and cupboard
void BuildingGen::generate()
{
  // create the building outer shape
  const double point_density = 600.0;
  const double building_width = random(7.0, 30.0);
  const double building_length = random(20.0, 80.0);
  const double building_height = random(3.0, 12.0);
  table_density = random(0.03, 0.07); // items per square metre
  cupboard_density = random(0.05, 0.15); // items per metre along wall

  Eigen::Vector3d building_pos(random(-20.0, 20.0), random(-20.0, 20.0), random(-20.0, 20.0));
  Cuboid building;
  building.min_bound = building_pos - 0.5*Eigen::Vector3d(building_width, building_length, building_height);
  building.max_bound = building_pos + 0.5*Eigen::Vector3d(building_width, building_length, building_height);

  // we generate the building as a list of negative cuboids, using a KD-tree type splitting
  std::vector<Cuboid> cuboids;
  std::vector<Cuboid> furniture;
  split(building, cuboids, furniture);

  // now we need to generate a path through the building...  
  std::vector<Eigen::Vector3d> points;
  std::vector<int> visited(cuboids.size());
  for (int i = 0; i<(int)visited.size(); i++)//auto &cuboid: cuboids)
    visited[i] = 0;
  buildPath((cuboids[0].max_bound + cuboids[0].min_bound)/2.0, 0, -1, cuboids, points, point_density, visited);
  std::vector<Eigen::Vector3d> dirs;
  for (int i = 0; i<(int)points.size(); i++)
    dirs.push_back(Eigen::Vector3d(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0)).normalized());

  // lets add a few windows and an outside area
  const double big = 1e4;
  const double wall_width = 0.3;
  const double window_width = 1.2;
  const double window_height = 1.5;
  Cuboid right = building; 
  right.min_bound[0] = building.max_bound[0] + wall_width;
  right.max_bound[0] = big;
  cuboids.push_back(right);
  Cuboid left = building; 
  left.max_bound[0] = building.min_bound[0] - wall_width;
  left.min_bound[0] = -big;
  cuboids.push_back(left);
  int num_windows = rand()%11;
  for (int i = 0; i<num_windows; i++)
  {
    Eigen::Vector3d pos = building.min_bound + Eigen::Vector3d(0, building_length*(double)(i+1)/(double)(num_windows+1), window_height);
    Cuboid window;
    window.min_bound = pos - 0.5*Eigen::Vector3d(window_width, window_width, window_width);
    window.max_bound = pos + 0.5*Eigen::Vector3d(window_width, window_width, window_width);
    cuboids.push_back(window);
    window.min_bound[0] += building_width;
    window.max_bound[0] += building_width;
    cuboids.push_back(window);
  }

  // now for every ray, we need to intersect it with the cuboids...
  // currently brute force, for p points and c cuboids it takes O(pc^2) because each cuboid has to test overlap with every other

  // alternative1: maintain a sorted list of in/out locations, for each cuboid, insert the start/end points and remove excess
  // should be O(pcn) where n are the number of actually intersecting cuboids, or O(pclog n) if a binary insert

  // alternative2: create a graph of overlapping negative spaces, walk the graph, (searching the start node to start with)
  // should be: O(p)  <-- this looks not much harder than alternative1, but definitely faster

  for (int i = 0; i<(int)points.size(); i++)
  {
    Eigen::Vector3d &start = points[i];
    Eigen::Vector3d &dir = dirs[i];

    const double max_range = 20.0;
    std::vector<Eigen::Vector3d> hits;
    for (int i = 0; i < (int)cuboids.size(); i++)
    {
      double new_range = max_range;
      if (cuboids[i].rayIntersectNegativeBox(start, dir, new_range))
        hits.push_back(start + dir * (new_range + 1e-6));
    }
    double range = max_range;
    for (auto &hit : hits)
    {
      bool intersected = false;
      for (auto &cuboid : cuboids)
      {
        if (cuboid.intersects(hit))
        {
          intersected = true;
          break;
        }
      }
      if (!intersected)
        range = std::min(range, (hit - start).norm());
    }
    for (auto &cuboid : furniture) 
      cuboid.rayIntersectBox(start, dir, range);

    const double range_noise = 0.03;
    Eigen::Vector3d ray_end = start + (range + random(-range_noise, range_noise)) * dir;
    ray_starts_.push_back(start);
    ray_ends_.push_back(ray_end);
    ray_bounded_.push_back(range != max_range);
  }
}
} // ray