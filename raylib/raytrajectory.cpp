// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrajectory.h"

namespace ray
{
void Trajectory::save(const std::string &file_name, double time_offset)
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  ofs.unsetf(std::ios::floatfield);
  ofs.precision(15);
  ofs << "%time x y z q0 q1 q2 q3 userfields" << std::endl;
  for (size_t i = 0; i < nodes.size(); i++)
  {
    const Pose &pose = nodes[i].pose;
    ofs << nodes[i].time + time_offset << " " << pose.position[0] << " " << pose.position[1] << " " << pose.position[2]
        << " " << pose.rotation.w() << " " << pose.rotation.x() << " " << pose.rotation.y() << " " << pose.rotation.z()
        << " " << std::endl;
  }
}

/**Loads the trajectory into the supplied vector and returns if successful*/
bool Trajectory::load(const std::string &file_name)
{
  std::cout << "loading " << file_name << std::endl;
  std::string line;
  int size = -1;
  {
    std::ifstream ifs(file_name.c_str(), std::ios::in);
    if (!ifs)
    {
      std::cerr << "Failed to open trajectory file: " << file_name << std::endl;
      return false;
    }
    ASSERT(ifs.is_open());
    getline(ifs, line);

    while (!ifs.eof())
    {
      getline(ifs, line);
      size++;
    }
  }
  std::ifstream ifs(file_name.c_str(), std::ios::in);
  if (!ifs)
  {
    std::cerr << "Failed to open trajectory file: " << file_name << std::endl;
    return false;
  }
  getline(ifs, line);
  std::vector<Node> trajectory(size);
  for (size_t i = 0; i < trajectory.size(); i++)
  {
    if (!ifs)
    {
      std::cerr << "Invalid stream when loading trajectory file: " << file_name << std::endl;
      return false;
    }

    getline(ifs, line);
    std::istringstream iss(line);
    iss >> trajectory[i].time >> trajectory[i].pose.position[0] >> trajectory[i].pose.position[1] >>
      trajectory[i].pose.position[2] >> trajectory[i].pose.rotation.w() >> trajectory[i].pose.rotation.x() >>
      trajectory[i].pose.rotation.y() >> trajectory[i].pose.rotation.z();
  }

  if (!ifs)
  {
    std::cerr << "Invalid stream when loading trajectory file: " << file_name << std::endl;
    return false;
  }

  // All is well add the loaded trajectory in
  nodes = trajectory;

  return true;
}
} // ray