// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raybuildinggen.h"

namespace ray
{
void BuildingGen::splitRoom(const Cuboid &cuboid, std::vector<Cuboid> &cuboids,
                            std::vector<std::vector<Cuboid>> &furniture)
{
  Eigen::Vector3d ext = cuboid.max_bound_ - cuboid.min_bound_;
  double sum = ext[0] + ext[1] + ext[2];
  double r = random(0.0, sum);
  int ax = r < ext[0] ? 0 : r < ext[0] + ext[1] ? 1 : 2;
  if (ext[2] > 2.0 * params_.room_scales[2] && random(0.0, 1.0) < params_.distinct_floor_likelihood)
    ax = 2.0;
  double length = ext[ax];
  if (length > 2.0 * params_.room_scales[ax])
  {
    double split_point = cuboid.min_bound_[ax] + random(params_.room_scales[ax], length - params_.room_scales[ax]);
    Cuboid rooms[2] = { cuboid, cuboid };
    rooms[0].max_bound_[ax] = split_point;
    rooms[1].min_bound_[ax] = split_point;

    if (ax == 2)  // stair hole
    {
      // just a hole for now
      Cuboid hole;
      Eigen::Vector3d hole_pos = cuboid.min_bound_ + ext * random(0.25, 0.75);
      hole_pos[ax] = split_point;
      hole.min_bound_ =
        hole_pos - Eigen::Vector3d(1.5 * params_.door_width, 2.5 * params_.door_width, params_.door_width);
      hole.max_bound_ =
        hole_pos + Eigen::Vector3d(1.5 * params_.door_width, 2.5 * params_.door_width, params_.door_width);
      cuboids.push_back(hole);
    }
    else
    {
      Cuboid door;
      Eigen::Vector3d door_pos = cuboid.min_bound_ + ext * random(0.25, 0.75);
      door_pos[ax] = split_point;
      door.min_bound_ = door_pos - Eigen::Vector3d(params_.door_width, params_.door_width, 0.0);
      door.max_bound_ = door_pos + Eigen::Vector3d(params_.door_width, params_.door_width, 0.0);
      door.min_bound_[2] = cuboid.min_bound_[2];
      door.max_bound_[2] = cuboid.min_bound_[2] + params_.door_height;
      cuboids.push_back(door);
    }
    furniture.push_back(std::vector<Cuboid>());
    splitRoom(rooms[0], cuboids, furniture);
    splitRoom(rooms[1], cuboids, furniture);
  }
  else
  {
    Cuboid room = cuboid;
    room.min_bound_ += 0.5 * Eigen::Vector3d(params_.wall_width, params_.wall_width, 0);
    room.max_bound_ -= 0.5 * Eigen::Vector3d(params_.wall_width, params_.wall_width, 2.0 * params_.floor_width);

    // now lets add some furniture...
    // tables:
    double floor_area = (room.max_bound_[0] - room.min_bound_[0]) * (room.max_bound_[1] - room.min_bound_[1]);
    int num_tables = (int)(floor_area * params_.table_density);
    std::vector<Cuboid> added;
    for (int i = 0; i < num_tables; i++)
    {
      double table_width = 1.0;
      double table_length = 1.5;
      const double table_height = 1.1;
      if (ray::rand() % 2)
        std::swap(table_width, table_length);
      Eigen::Vector3d table_rad = 0.5 * Eigen::Vector3d(table_width, table_length, table_height);
      Eigen::Vector3d table_pos(
        random(room.min_bound_[0] + table_width / 2.0, room.max_bound_[0] - table_width / 2.0),
        random(room.min_bound_[1] + table_length / 2.0, room.max_bound_[1] - table_length / 2.0),
        room.min_bound_[2] + table_height / 2.0);
      Cuboid table;
      table.min_bound_ = table_pos - table_rad;
      table.max_bound_ = table_pos + table_rad;
      added.push_back(table);
    }
    // cupboards:
    double wall_length = 2.0 * ((room.max_bound_[0] - room.min_bound_[0]) + (room.max_bound_[1] - room.min_bound_[1]));
    int num_cupboards = (int)(wall_length * params_.cupboard_density);
    for (int i = 0; i < num_cupboards; i++)
    {
      double cupboard_width = 0.6;
      double cupboard_length = 2.0;
      double cupboard_height = random(1.0, 3.0);
      int wall_id = ray::rand() % 2;
      int side = ray::rand() % 2;
      if (wall_id == 1)
        std::swap(cupboard_width, cupboard_length);
      Eigen::Vector3d cupboard_rad = 0.5 * Eigen::Vector3d(cupboard_width, cupboard_length, cupboard_height);
      Eigen::Vector3d cupboard_pos;
      cupboard_pos[2] = room.min_bound_[2] + 0.5 * cupboard_height;
      if (wall_id == 0)
      {
        cupboard_pos[0] = side ? room.min_bound_[0] + 0.5 * cupboard_width : room.max_bound_[0] - 0.5 * cupboard_width;
        cupboard_pos[1] =
          random(room.min_bound_[1] + 0.5 * cupboard_length, room.max_bound_[1] - 0.5 * cupboard_length);
      }
      else
      {
        cupboard_pos[0] = random(room.min_bound_[0] + 0.5 * cupboard_width, room.max_bound_[0] - 0.5 * cupboard_width);
        cupboard_pos[1] =
          side ? room.min_bound_[1] + 0.5 * cupboard_length : room.max_bound_[1] - 0.5 * cupboard_length;
      }
      Cuboid cupboard;
      cupboard.min_bound_ = cupboard_pos - cupboard_rad;
      cupboard.max_bound_ = cupboard_pos + cupboard_rad;
      added.push_back(cupboard);
    }
    bool overlaps = false;
    size_t f_id = cuboids.size();
    furniture.push_back(std::vector<Cuboid>());
    for (auto &item : added)
    {
      for (auto &other : furniture[f_id])
        if (other.overlaps(item))
          overlaps = true;
      for (auto &other : cuboids)
        if (other.overlaps(item))
          overlaps = true;
      if (!overlaps)
        furniture[f_id].push_back(item);
    }
    cuboids.push_back(room);
  }
}

void buildPath(const Eigen::Vector3d &start_pos, size_t id, size_t last_id, const std::vector<Cuboid> &cuboids,
               std::vector<Eigen::Vector3d> &points, double density, std::vector<int> &visited)
{
  visited[id] = 1;
  Eigen::Vector3d last_pos = start_pos;
  for (size_t i = 0; i < cuboids.size(); i++)
  {
    if (i == id || i == last_id || visited[i])
      continue;
    if (cuboids[i].overlaps(cuboids[id]))
    {
      Eigen::Vector3d pos = (cuboids[i].max_bound_ + cuboids[i].min_bound_) / 2.0;

      int num = (int)((pos - last_pos).norm() * density);
      for (int j = 0; j < num; j++) points.push_back(last_pos + (pos - last_pos) * (double)j / (double)num);
      buildPath(pos, i, id, cuboids, points, density, visited);
      last_pos = pos;
    }
  }
  if (last_pos == start_pos)  // no outlets found
  {
    Eigen::Vector3d pos = (cuboids[id].max_bound_ + cuboids[id].min_bound_) / 2.0;
    int num = (int)((pos - last_pos).norm() * density);
    for (int j = 0; j < num; j++) points.push_back(last_pos + (pos - last_pos) * (double)j / (double)num);
    last_pos = pos;
  }
  int num = (int)((start_pos - last_pos).norm() * density);
  for (int j = 0; j < num; j++) points.push_back(last_pos + (start_pos - last_pos) * (double)j / (double)num);
}

// generates a set of rooms within a cuboidal frame
void BuildingGen::generate()
{
  std::cout << "generating" << std::endl;
  // create the building outer shape
  const double point_density = 3000.0;
  const double building_width = random(7.0, 30.0);
  const double building_length = random(20.0, 80.0);
  const double building_height = random(3.0, 14.0);
  params_.table_density = random(0.03, 0.07);     // items per square metre
  params_.cupboard_density = random(0.05, 0.15);  // items per metre along wall

  Eigen::Vector3d building_pos(random(-20.0, 20.0), random(-20.0, 20.0), random(-20.0, 20.0));
  Cuboid building;
  building.min_bound_ = building_pos - 0.5 * Eigen::Vector3d(building_width, building_length, building_height);
  building.max_bound_ = building_pos + 0.5 * Eigen::Vector3d(building_width, building_length, building_height);

  // we generate the building as a list of negative cuboids, using a KD-tree type splitting
  std::vector<Cuboid> cuboids;
  std::vector<std::vector<Cuboid>> furniture;  // each room has furnitude added (also cuboids)
  splitRoom(building, cuboids, furniture);

  // now we need to generate a path through the building...
  std::vector<Eigen::Vector3d> points;
  std::vector<int> visited(cuboids.size());  // keep track of rooms we have visited on the path through
  for (size_t i = 0; i < visited.size(); i++) visited[i] = 0;
  buildPath((cuboids[0].max_bound_ + cuboids[0].min_bound_) / 2.0, 0, -1, cuboids, points, point_density, visited);
  // give each path point a random ray direction
  std::vector<Eigen::Vector3d> dirs;
  for (size_t i = 0; i < points.size(); i++)
    dirs.push_back(Eigen::Vector3d(random(-1.0, 1.0), random(-1.0, 1.0), random(-1.0, 1.0)).normalized());

  // add a few windows and an outside area
  const double big = 1e4;
  Cuboid right = building;
  right.min_bound_[0] = building.max_bound_[0] + params_.outer_wall_width;
  right.max_bound_[0] = big;
  cuboids.push_back(right);
  Cuboid left = building;
  left.max_bound_[0] = building.min_bound_[0] - params_.outer_wall_width;
  left.min_bound_[0] = -big;
  cuboids.push_back(left);
  int num_windows = ray::rand() % 11;
  for (int i = 0; i < num_windows; i++)
  {
    Eigen::Vector3d pos =
      building.min_bound_ +
      Eigen::Vector3d(0, building_length * (double)(i + 1) / (double)(num_windows + 1), params_.window_height);
    Cuboid window;
    window.min_bound_ = pos - 0.5 * Eigen::Vector3d(params_.window_width, params_.window_width, params_.window_width);
    window.max_bound_ = pos + 0.5 * Eigen::Vector3d(params_.window_width, params_.window_width, params_.window_width);
    cuboids.push_back(window);
    window.min_bound_[0] += building_width;
    window.max_bound_[0] += building_width;
    cuboids.push_back(window);
  }
  furniture.resize(cuboids.size());

  std::cout << "ray casting" << std::endl;

  // create a graph of overlapping negative spaces, walk the graph, (searching the start node to start with)
  // should be: O(p)  <-- this looks not much harder than alternative1, but definitely faster
  std::vector<std::vector<size_t>> neighbours(cuboids.size());
  for (size_t i = 0; i < cuboids.size(); i++)
    for (size_t j = 0; j < cuboids.size(); j++)
      if (i != j && cuboids[i].overlaps(cuboids[j]))
        neighbours[i].push_back(j);

  size_t start_id = 0;  // if this is not the right initial cuboid, it will search to find the right one
  // for each ray (points, dirs) intersect it with the negative spaces (cuboids) and the furniture to get a range
  for (size_t i = 0; i < points.size(); i++)
  {
    Eigen::Vector3d &start = points[i];
    Eigen::Vector3d &dir = dirs[i];

    const double max_range = 20.0;
    double range = max_range;
    // first, adjust startID:
    if (!cuboids[start_id].intersects(start))  // if it doesn't intersect then search for an intersecting room
    {
      bool found_start = false;
      for (auto &id : neighbours[start_id])  // initially try the neighbours
      {
        if (cuboids[id].intersects(start))
        {
          start_id = id;
          found_start = true;
          break;
        }
      }
      if (!found_start)  // otherwise a brute force search
      {
        for (size_t id = 0; id < cuboids.size(); id++)
        {
          if (cuboids[id].intersects(start))
          {
            start_id = id;
            found_start = true;
            break;
          }
        }
        if (!found_start)  // failing that, the ray must start outside the building, so don't add it
          continue;
      }
    }
    // get range in box it already intersects
    cuboids[start_id].intersectsRay(start, dir, range, false);
    const double eps = 1e-6;
    Eigen::Vector3d hit = start + dir * (range + eps);

    size_t id = start_id;
    bool found = range < max_range;
    for (auto &cuboid : furniture[id])  // if the room furniture is closer, then use this
    {
      if (cuboid.intersectsRay(start, dir, range, true))
      {
        found = false;
        break;
      }
    }
    // now check adjacent rooms (and doorways, windows) to see if the ray can travel further
    while (found)  // Note, this can't cause an infinite cycle because rays are straight lines
    {
      found = false;
      for (auto &other_id : neighbours[id])
      {
        if (!cuboids[other_id].intersects(hit))
          continue;
        // hit intersects, so we extend the ray as far as this cuboid reaches
        double new_range = max_range;
        cuboids[other_id].intersectsRay(start, dir, new_range, false);
        range = new_range;
        hit = start + dir * (range + eps);
        id = other_id;
        found = new_range < max_range;
        for (auto &cuboid : furniture[id])  // also check the furniture in this space
        {
          if (cuboid.intersectsRay(start, dir, range, true))
          {
            found = false;
            break;
          }
        }
        break;
      }
    }

    // finally we have a range value for the ray, so add it to the ray cloud
    const double range_noise = 0.03;
    Eigen::Vector3d ray_end = start + (range + random(-range_noise, range_noise)) * dir;
    ray_starts_.push_back(start);
    ray_ends_.push_back(ray_end);
    ray_bounded_.push_back(range != max_range);
  }
}
}  // namespace ray