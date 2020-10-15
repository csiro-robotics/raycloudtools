// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrajectory.h"

namespace ray
{
void Trajectory::save(const std::string &file_name)
{
  std::cout << "saving trajectory " << file_name << std::endl;
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  ofs.unsetf(std::ios::floatfield);
  ofs.precision(15);
  ofs << "%time x y z userfields" << std::endl;
  for (size_t i = 0; i < nodes.size(); i++)
  {
    const Eigen::Vector3d &pos = nodes[i].pos;
    ofs << nodes[i].time << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << std::endl;
  }
}

/**Loads the trajectory into the supplied vector and returns if successful*/
bool Trajectory::load(const std::string &file_name)
{
  std::cout << "loading trajectory " << file_name << std::endl;
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
    iss >> trajectory[i].time >> trajectory[i].pos[0] >> trajectory[i].pos[1] >>
      trajectory[i].pos[2];
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